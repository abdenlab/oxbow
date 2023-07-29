use std::fs::File;
use std::io::{self, Read, Seek};
use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, Float32Builder, GenericStringBuilder, Int32Builder,
    StringArray, StringDictionaryBuilder,
};
use arrow::{
    datatypes::Int32Type, error::ArrowError, record_batch::RecordBatch,
};
use noodles::core::Region;
use noodles::{bgzf, csi, tabix, vcf};

use crate::batch_builder::{write_ipc, BatchBuilder};

type BufferedReader = io::BufReader<File>;


/// A VCF reader.
pub struct VcfReader {
    reader: vcf::Reader<bgzf::Reader<BufferedReader>>,
    header: vcf::Header,
    index: csi::Index,
}

impl VcfReader {
    /// Creates a VCF Reader.
    pub fn new(path: &str) -> std::io::Result<Self> {
        let tbi_path = format!("{}.tbi", path);
        let csi_path = format!("{}.csi", path);
        let index = if Path::new(&tbi_path).exists() {
            tabix::read(tbi_path)?
        } else if Path::new(&csi_path).exists() {
            csi::read(csi_path)?
        } else {
            panic!("Could not find a .tbi or .csi index file for the given VCF file.");
        };

        let file = std::fs::File::open(path)?;
        let buf_file = std::io::BufReader::with_capacity(1024 * 1024, file);
        let mut reader = vcf::Reader::new(bgzf::Reader::new(buf_file));
        let header = reader.read_header()?;
        Ok(Self { reader, header, index })
    }

    /// Returns the records in the given region as Apache Arrow IPC.
    ///
    /// If the region is `None`, all records are returned.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use oxbow::vcf::VcfReader;
    ///
    /// let mut reader = VcfReader::new("sample.vcf.gz").unwrap();
    /// let ipc = reader.records_to_ipc(Some("sq0:1-1000")).unwrap();
    /// ```
    pub fn records_to_ipc(&mut self, region: Option<&str>) -> Result<Vec<u8>, ArrowError> {
        let batch_builder = VcfBatchBuilder::new(1024, &self.header)?;
        if let Some(region) = region {
            let region: Region = region.parse().unwrap();
            let query = self
                .reader
                .query(&self.header, &self.index, &region)
                .unwrap()
                .map(|r| r.unwrap());
            return write_ipc(query, batch_builder);
        }
        let records = self.reader.records(&self.header).map(|r| r.unwrap());
        write_ipc(records, batch_builder)
    }

    pub fn records_to_ipc_from_vpos(&mut self, pos_lo: (u64, u16), pos_hi: (u64, u16)) -> Result<Vec<u8>, ArrowError> {
        let vpos_lo = bgzf::VirtualPosition::try_from(pos_lo).unwrap();
        let vpos_hi = bgzf::VirtualPosition::try_from(pos_hi).unwrap();
        let batch_builder = VcfBatchBuilder::new(1024, &self.header)?;
        let records = VcfRecords::new(&mut self.reader, &self.header, vpos_lo, vpos_hi).map(|r| r.unwrap());
        write_ipc(records, batch_builder)
    }
}

struct VcfBatchBuilder {
    chrom: StringDictionaryBuilder<Int32Type>,
    pos: Int32Builder,
    id: GenericStringBuilder<i32>,
    ref_: GenericStringBuilder<i32>,
    alt: GenericStringBuilder<i32>,
    qual: Float32Builder,
    filter: GenericStringBuilder<i32>,
    info: GenericStringBuilder<i32>,
    format: GenericStringBuilder<i32>,
}

impl VcfBatchBuilder {
    pub fn new(capacity: usize, header: &vcf::Header) -> Result<Self, ArrowError> {
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
        })
    }
}

impl BatchBuilder for VcfBatchBuilder {
    type Record = vcf::record::Record;

    fn push(&mut self, record: &Self::Record) {
        self.chrom.append_value(record.chromosome().to_string());
        self.pos.append_value(usize::from(record.position()) as i32);
        self.id.append_value(record.ids().to_string());
        self.ref_.append_value(record.reference_bases().to_string());
        self.alt.append_value(record.alternate_bases().to_string());
        self.qual
            .append_option(record.quality_score().map(f32::from));
        self.filter
            .append_option(record.filters().map(|f| f.to_string()));
        self.info.append_value(record.info().to_string());
        self.format.append_value(record.format().to_string());
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


// Reads VCF Records from a virtualposition range in a VCF file
pub struct VcfRecords<'a, R> {
    reader: &'a mut vcf::Reader<bgzf::reader::Reader<R>>,
    header: &'a vcf::Header,
    record: vcf::record::Record,
    vpos_lo: bgzf::VirtualPosition,
    vpos_hi: bgzf::VirtualPosition,
}

impl<'a, R> VcfRecords<'a, R>
where
    R: Read + Seek,
{
    pub(crate) fn new(
        reader: &'a mut vcf::Reader<bgzf::reader::Reader<R>>, 
        header: &'a vcf::Header,
        vpos_lo: bgzf::VirtualPosition,
        vpos_hi: bgzf::VirtualPosition,
    ) -> Self {
        let _ = reader.seek(vpos_lo);
        Self {
            reader,
            header,
            record: vcf::Record::default(),
            vpos_lo,
            vpos_hi,
        }
    }

    pub fn reset(&mut self) -> Option<bgzf::VirtualPosition> {
        self.reader.seek(self.vpos_lo).ok()
    }
}

impl<'a, R> Iterator for VcfRecords<'a, R>
where
    R: Read + Seek,
{
    type Item = io::Result<vcf::Record>;
    
    fn next(&mut self) -> Option<Self::Item> {
        if self.reader.virtual_position() >= self.vpos_hi {
            return None;
        }    

        match self.reader.read_record(self.header, &mut self.record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.record.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}
