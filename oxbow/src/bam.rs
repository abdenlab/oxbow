use std::fs::File;
use std::io::{self, Read, Seek};
use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, GenericStringBuilder, Int32Array, Int32Builder, StringArray, StringDictionaryBuilder,
    UInt16Array, UInt16Builder, UInt8Array, UInt8Builder,
};
use arrow::{
    datatypes::Int32Type, error::ArrowError, record_batch::RecordBatch,
};
use noodles::core::Region;
use noodles::{bam, bgzf, csi, sam};

use crate::batch_builder::{write_ipc, BatchBuilder};

type BufferedReader = io::BufReader<File>;


/// A BAM reader.
pub struct BamReader {
    reader: bam::Reader<bgzf::Reader<BufferedReader>>,
    header: sam::Header,
    index: csi::Index,
}

impl BamReader {
    /// Creates a BAM reader.
    pub fn new(path: &str) -> std::io::Result<Self> {
        let bai_path = format!("{}.bai", path);
        let csi_path = format!("{}.csi", path);
        let index = if Path::new(&bai_path).exists() {
            bam::bai::read(bai_path)?
        } else if Path::new(&csi_path).exists() {
            csi::read(csi_path)?
        } else {
            panic!("Could not find a .bai or .csi index file for the given BAM file.");
        };

        let file = std::fs::File::open(path)?;
        let buf_file = std::io::BufReader::with_capacity(1024 * 1024, file);
        let mut reader = bam::Reader::new(buf_file);
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
    /// use oxbow::bam::BamReader;
    ///
    /// let mut reader = BamReader::new("sample.bam").unwrap();
    /// let ipc = reader.records_to_ipc(Some("sq0:1-1000")).unwrap();
    /// ```
    pub fn records_to_ipc(&mut self, region: Option<&str>) -> Result<Vec<u8>, ArrowError> {
        let batch_builder = BamBatchBuilder::new(1024, &self.header)?;
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
        let batch_builder = BamBatchBuilder::new(1024, &self.header)?;
        let records = BamRecords::new(&mut self.reader, &self.header, vpos_lo, vpos_hi).map(|r| r.unwrap());
        write_ipc(records, batch_builder)
    }
}

struct BamBatchBuilder<'a> {
    header: &'a sam::Header,
    qname: GenericStringBuilder<i32>,
    flag: UInt16Builder,
    rname: StringDictionaryBuilder<Int32Type>,
    pos: Int32Builder,
    mapq: UInt8Builder,
    cigar: GenericStringBuilder<i32>,
    rnext: StringDictionaryBuilder<Int32Type>,
    pnext: Int32Builder,
    tlen: Int32Builder,
    seq: GenericStringBuilder<i32>,
    qual: GenericStringBuilder<i32>,
    end: Int32Builder,
}

impl<'a> BamBatchBuilder<'a> {
    pub fn new(capacity: usize, header: &'a sam::Header) -> Result<Self, ArrowError> {
        let categories = StringArray::from(
            header
                .reference_sequences()
                .iter()
                .map(|(rs, _)| Some(rs.as_str()))
                .collect::<Vec<_>>(),
        );
        Ok(Self {
            header,
            qname: GenericStringBuilder::<i32>::new(),
            flag: UInt16Array::builder(capacity),
            rname: StringDictionaryBuilder::<Int32Type>::new_with_dictionary(
                capacity,
                &categories,
            )?,
            pos: Int32Array::builder(capacity),
            mapq: UInt8Array::builder(capacity),
            cigar: GenericStringBuilder::<i32>::new(),
            rnext: StringDictionaryBuilder::<Int32Type>::new_with_dictionary(
                capacity,
                &categories,
            )?,
            pnext: Int32Array::builder(capacity),
            tlen: Int32Array::builder(capacity),
            seq: GenericStringBuilder::<i32>::new(),
            qual: GenericStringBuilder::<i32>::new(),
            end: Int32Array::builder(capacity),
        })
    }
}

impl<'a> BatchBuilder for BamBatchBuilder<'a> {
    type Record = sam::alignment::Record;

    fn push(&mut self, record: &Self::Record) {
        self.qname.append_option(record.read_name());
        self.flag.append_value(record.flags().bits());
        let rname = match record.reference_sequence(self.header) {
            Some(Ok((name, _))) => Some(name.as_str()),
            _ => None,
        };
        self.rname.append_option(rname);
        self.pos
            .append_option(record.alignment_start().map(|x| x.get() as i32));
        self.mapq
            .append_option(record.mapping_quality().map(|x| x.get()));
        self.cigar.append_value(record.cigar().to_string());
        let rnext = match record.mate_reference_sequence(self.header) {
            Some(Ok((name, _))) => Some(name.as_str()),
            _ => None,
        };
        self.rnext.append_option(rnext);
        self.pnext
            .append_option(record.mate_alignment_start().map(|x| x.get() as i32));
        self.tlen.append_value(record.template_length());
        self.seq.append_value(record.sequence().to_string());
        self.qual.append_value(record.quality_scores().to_string());

        // extra
        self.end
            .append_option(record.alignment_end().map(|x| x.get() as i32));
    }

    fn finish(mut self) -> Result<RecordBatch, ArrowError> {
        RecordBatch::try_from_iter(vec![
            // spec
            ("qname", Arc::new(self.qname.finish()) as ArrayRef),
            ("flag", Arc::new(self.flag.finish()) as ArrayRef),
            ("rname", Arc::new(self.rname.finish()) as ArrayRef),
            ("pos", Arc::new(self.pos.finish()) as ArrayRef),
            ("mapq", Arc::new(self.mapq.finish()) as ArrayRef),
            ("cigar", Arc::new(self.cigar.finish()) as ArrayRef),
            ("rnext", Arc::new(self.rnext.finish()) as ArrayRef),
            ("pnext", Arc::new(self.pnext.finish()) as ArrayRef),
            ("tlen", Arc::new(self.tlen.finish()) as ArrayRef),
            ("seq", Arc::new(self.seq.finish()) as ArrayRef),
            ("qual", Arc::new(self.qual.finish()) as ArrayRef),
            // extra
            ("end", Arc::new(self.end.finish()) as ArrayRef),
        ])
    }
}


// Reads SAM records from a virtualposition range in a BAM file
pub struct BamRecords<'a, R>
where
    R: Read + Seek,
{
    reader: &'a mut bam::Reader<bgzf::reader::Reader<R>>,
    header: &'a sam::Header,
    record: sam::alignment::Record,
    vpos_lo: bgzf::VirtualPosition,
    vpos_hi: bgzf::VirtualPosition,
}

impl<'a, R> BamRecords<'a, R>
where
    R: Read + Seek,
{
    pub fn new(
        reader: &'a mut bam::Reader<bgzf::reader::Reader<R>>,
        header: &'a sam::Header, 
        vpos_lo: bgzf::VirtualPosition, 
        vpos_hi: bgzf::VirtualPosition
    ) -> Self {
        let _ = reader.seek(vpos_lo);
        Self {
            reader,
            header,
            record: sam::alignment::Record::default(),
            vpos_lo,
            vpos_hi,
        }
    }

    pub fn reset(&mut self) -> Option<bgzf::VirtualPosition> {
        self.reader.seek(self.vpos_lo).ok()
    }
}

impl<'a, R> Iterator for BamRecords<'a, R>
where
    R: Read + Seek,
{
    type Item = io::Result<sam::alignment::Record>;

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


#[cfg(test)]
mod tests {
    use super::*;
    use arrow::ipc::reader::FileReader;
    use arrow::record_batch::RecordBatch;

    fn read_record_batch(region: Option<&str>) -> RecordBatch {
        let mut dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push("../fixtures/sample.bam");
        let mut reader = BamReader::new(dir.to_str().unwrap()).unwrap();
        let ipc = reader.records_to_ipc(region).unwrap();
        let cursor = std::io::Cursor::new(ipc);
        let mut arrow_reader = FileReader::try_new(cursor, None).unwrap();
        // make sure we have one batch
        assert_eq!(arrow_reader.num_batches(), 1);
        arrow_reader.next().unwrap().unwrap()
    }

    #[test]
    fn test_read_all() {
        let record_batch = read_record_batch(None);
        assert_eq!(record_batch.num_rows(), 6);
    }

    #[test]
    fn test_region_full() {
        let record_batch = read_record_batch(Some("chr1"));
        assert_eq!(record_batch.num_rows(), 4);
    }

    #[test]
    fn rest_region_partial() {
        let record_batch = read_record_batch(Some("chr1:1-100000"));
        assert_eq!(record_batch.num_rows(), 2);
    }
}
