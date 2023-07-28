use arrow::array::{
    ArrayRef, GenericStringBuilder, Int32Array, Int32Builder, StringArray, StringDictionaryBuilder,
    UInt16Array, UInt16Builder, UInt8Array, UInt8Builder,
};
use arrow::{
    datatypes::Int32Type, error::ArrowError, record_batch::RecordBatch,
};
use noodles::core::Region;
use noodles::{bam, bgzf, sam};
use std::sync::Arc;

use crate::batch_builder::{write_ipc, BatchBuilder};
use crate::vpos;

type BufferedReader = std::io::BufReader<std::fs::File>;


/// A BAM reader.
pub struct BamReader {
    reader: bam::IndexedReader<bgzf::Reader<BufferedReader>>,
    unindexed_reader: bam::Reader<bgzf::Reader<std::fs::File>>,
    header: sam::Header,
}

impl BamReader {
    /// Creates a BAM reader.
    pub fn new(path: &str) -> std::io::Result<Self> {
        let index = bam::bai::read(format!("{}.bai", path))?;
        let file = std::fs::File::open(path)?;
        let bufreader = std::io::BufReader::with_capacity(1024 * 1024, file);
        let mut reader = bam::indexed_reader::Builder::default()
            .set_index(index)
            .build_from_reader(bufreader)?;
        let header = reader.read_header()?;
        let unindexed_reader = bam::reader::Builder::default().build_from_path(path).unwrap();
        Ok(Self { reader, unindexed_reader, header })
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
                .query(&self.header, &region)
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
        let records = vpos::BamRecords::new(&mut self.unindexed_reader, &self.header, vpos_lo, vpos_hi).map(|r| r.unwrap());
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
