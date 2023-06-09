use arrow::array::{
    ArrayRef, Float32Builder, GenericStringBuilder, Int32Builder,
    StringArray, StringDictionaryBuilder,
};
use arrow::{
    datatypes::Int32Type, error::ArrowError, record_batch::RecordBatch,
};
use byteorder::{LittleEndian, ReadBytesExt};
use noodles::core::Region;
use noodles::bcf::header::StringMaps;
use noodles::{bcf, bgzf, csi, vcf};
use std::sync::Arc;
use std::{ffi::CStr, io};
use std::io::{Read};

use crate::batch_builder::{write_ipc, BatchBuilder};

type BufferedReader = std::io::BufReader<std::fs::File>;


fn read_u8<R>(reader: &mut R) -> io::Result<u8>
where
    R: Read,
{
    let mut buf = [0; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

fn read_raw_header<R>(reader: &mut R) -> io::Result<String>
where
    R: Read,
{
    let l_text = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut buf = vec![0; l_text];
    reader.read_exact(&mut buf)?;

    CStr::from_bytes_with_nul(&buf)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .and_then(|c_header| {
            c_header
                .to_str()
                .map(|s| s.into())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
}

fn read_bcf_header<R>(reader: &mut R) -> io::Result<(vcf::Header, bcf::header::StringMaps)>
where
    R: Read,
{
    // read the MAGIC number
    let mut buf = [0; 3];
    reader.read_exact(&mut buf)?;

    let mut reader = reader;

    // read the version
    read_u8(&mut reader)?;
    read_u8(&mut reader)?;

    // read the raw header as a string
    let s = read_raw_header(&mut reader)?;

    // parse the VCF header and the string maps
    let header: vcf::Header = s.parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let string_maps: bcf::header::StringMaps = s.parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok((header, string_maps))
}


/// A BCF reader.
pub struct BcfReader {
    reader: bcf::Reader<bgzf::Reader<BufferedReader>>,
    header: vcf::Header,
    string_maps: StringMaps,
    index: csi::Index,
}

impl BcfReader {
    /// Creates a BCF Reader.
    pub fn new(path: &str) -> std::io::Result<Self> {
        let file = std::fs::File::open(path)?;
        let index = csi::read(format!("{}.csi", path))?;
        let bufreader = std::io::BufReader::with_capacity(1024 * 1024, file);
        let reader = bcf::Reader::new(bufreader);

        let mut bgzf_reader = std::fs::File::open(path).map(bgzf::Reader::new)?;
        let (header, string_maps) = read_bcf_header(&mut bgzf_reader)?;

        Ok(Self { reader, header, string_maps, index })
    }

    /// Returns the records in the given region as Apache Arrow IPC.
    ///
    /// If the region is `None`, all records are returned.
    ///
    /// # Examples
    ///
    /// ``no_run
    /// use oxbow::bcf::BcfReader;
    ///
    /// let mut reader = BcfReader::new("sample.bcf").unwrap();
    /// let ipc = reader.records_to_ipc(Some("sq0:1-1000")).unwrap();
    /// ```
    pub fn records_to_ipc(&mut self, region: Option<&str>) -> Result<Vec<u8>, ArrowError> {
        let s = self.string_maps.clone();
        let batch_builder = BcfBatchBuilder::new(1024, &self.header, &s)?;
        let string_maps = StringMaps::from(&self.header);
        if let Some(region) = region {
            let region: Region = region.parse().unwrap();
            let query = self
                .reader
                .query(string_maps.contigs(), &self.index, &region)
                .unwrap()
                .map(|r| r.unwrap());
            return write_ipc(query, batch_builder);
        }
        let records = self.reader.records().map(|r| r.unwrap());
        write_ipc(records, batch_builder)
    }
}

struct BcfBatchBuilder {
    chrom: StringDictionaryBuilder<Int32Type>,
    pos: Int32Builder,
    id: GenericStringBuilder<i32>,
    ref_: GenericStringBuilder<i32>,
    alt: GenericStringBuilder<i32>,
    qual: Float32Builder,
    filter: GenericStringBuilder<i32>,
    info: GenericStringBuilder<i32>,
    format: GenericStringBuilder<i32>,
    header: vcf::Header,
    string_maps: StringMaps,
}

impl BcfBatchBuilder {
    pub fn new(capacity: usize, header: &vcf::Header, string_maps: &StringMaps) -> Result<Self, ArrowError> {
        let categories = StringArray::from(
            header
                .contigs()
                .keys()
                .map(|k| k.to_string())
                .collect::<Vec<_>>(),
        );
        Ok(Self {
            chrom: StringDictionaryBuilder::<Int32Type>::new_with_dictionary(
                capacity,
                &categories,
            )?,
            pos: Int32Builder::with_capacity(capacity),
            id: GenericStringBuilder::<i32>::new(),
            ref_: GenericStringBuilder::<i32>::new(),
            alt: GenericStringBuilder::<i32>::new(),
            qual: Float32Builder::with_capacity(capacity),
            filter: GenericStringBuilder::<i32>::new(),
            info: GenericStringBuilder::<i32>::new(),
            format: GenericStringBuilder::<i32>::new(),
            header: header.clone(),
            string_maps: string_maps.clone(),
        })
    }
}

impl BatchBuilder for BcfBatchBuilder {
    type Record = bcf::record::Record;

    fn push(&mut self, record: &Self::Record) {

        let vcf_record = record.try_into_vcf_record(
            &self.header, &self.string_maps
        ).unwrap();

        self.chrom.append_value(vcf_record.chromosome().to_string());
        self.pos.append_value(usize::from(vcf_record.position()) as i32);
        self.id.append_value(vcf_record.ids().to_string());
        self.ref_.append_value(vcf_record.reference_bases().to_string());
        self.alt.append_value(vcf_record.alternate_bases().to_string());
        self.qual
            .append_option(vcf_record.quality_score().map(f32::from));
        self.filter
            .append_option(vcf_record.filters().map(|f| f.to_string()));
        self.info.append_value(vcf_record.info().to_string());
        self.format.append_value(vcf_record.format().to_string());
    }

    fn finish(mut self) -> Result<RecordBatch, ArrowError> {
        RecordBatch::try_from_iter(vec![
            // spec
            ("chrom", Arc::new(self.chrom.finish()) as ArrayRef),
            ("pos", Arc::new(self.pos.finish()) as ArrayRef),
            ("id", Arc::new(self.id.finish()) as ArrayRef),
            ("ref", Arc::new(self.ref_.finish()) as ArrayRef),
            ("alt", Arc::new(self.alt.finish()) as ArrayRef),
            ("qual", Arc::new(self.qual.finish()) as ArrayRef),
            ("filter", Arc::new(self.filter.finish()) as ArrayRef),
            ("info", Arc::new(self.info.finish()) as ArrayRef),
            ("format", Arc::new(self.format.finish()) as ArrayRef),
        ])
    }
}
