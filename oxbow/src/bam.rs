use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, Read, Seek};
use std::path::Path;
use std::sync::Arc;

use arrow::array::ArrayBuilder;
use arrow::array::{
    ArrayRef, GenericStringBuilder, Int32Array, Int32Builder, StringArray, StringDictionaryBuilder,
    UInt16Array, UInt16Builder, UInt8Array, UInt8Builder,
};
use arrow::{datatypes::Int32Type, error::ArrowError, record_batch::RecordBatch};
use noodles::core::Region;
use noodles::sam::record::data::field::Tag;
use noodles::{bam, bgzf, csi, sam};
use std::collections::HashSet;
use std::str::FromStr;

use crate::batch_builder::{write_ipc_err, BatchBuilder};

pub fn index_from_reader<R>(mut read: R) -> io::Result<csi::Index>
where
    R: Read + Seek,
{
    // Unlike .tbi and .csi, .bai is not bgzf-compressed
    // so we read off the magic directly.
    let mut magic = [0; 4];
    read.read_exact(&mut magic)?;
    read.seek(io::SeekFrom::Start(0))?;
    if magic == b"BAI\x01" as &[u8] {
        let mut bai_reader = bam::bai::Reader::new(read);
        bai_reader.read_header()?;
        bai_reader.read_index()
    } else {
        let mut csi_reader = csi::Reader::new(read);
        csi_reader.read_index()
    }
}

pub fn index_from_path(path: &str) -> io::Result<csi::Index> {
    let bai_path = format!("{}.bai", path);
    let csi_path = format!("{}.csi", path);
    let index = if Path::new(&bai_path).exists() {
        bam::bai::read(bai_path)?
    } else if Path::new(&csi_path).exists() {
        csi::read(csi_path)?
    } else {
        panic!("Could not find a .bai or .csi index file for the given BAM file.");
    };
    Ok(index)
}

/// A BAM reader.
pub struct BamReader<R> {
    reader: bam::Reader<bgzf::Reader<R>>,
    header: sam::Header,
    index: csi::Index,
}

impl BamReader<BufReader<File>> {
    /// Creates a BAM reader from a given file path.
    pub fn new_from_path(path: &str) -> std::io::Result<Self> {
        let index = index_from_path(path)?;
        let file = std::fs::File::open(path)?;
        let buf_file = std::io::BufReader::with_capacity(1024 * 1024, file);
        let mut reader = bam::Reader::new(buf_file);
        let header = reader.read_header()?;

        Ok(Self {
            reader,
            header,
            index,
        })
    }
}

impl<R: Read + Seek> BamReader<R> {
    /// Creates a BAM reader.
    pub fn new(read: R, index: csi::Index) -> std::io::Result<Self> {
        let mut reader = bam::Reader::new(read);
        let header = reader.read_header()?;

        Ok(Self {
            reader,
            header,
            index,
        })
    }

    /// Returns the records in the given region as Apache Arrow IPC.
    ///
    /// If the region is `None`, all records are returned. The second paramter to 
    /// `records_to_ipc` is a set of tags to include in the output. If it is `None`,
    /// all tags are included. If it is `Some`, only the tags in the set are included.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use oxbow::bam::BamReader;
    ///
    /// let mut reader = BamReader::new_from_path("sample.bam").unwrap();
    /// let ipc = reader.records_to_ipc(Some("sq0:1-1000"), None).unwrap();
    /// ```
    pub fn records_to_ipc(
        &mut self,
        region: Option<&str>,
        tags: Option<HashSet<&str>>,
    ) -> Result<Vec<u8>, ArrowError> {
        let batch_builder = BamBatchBuilder::new(1024, &self.header, tags)?;
        if let Some(region) = region {
            let region: Region = region.parse().unwrap();
            let query = self
                .reader
                .query(&self.header, &self.index, &region)
                .map_err(|e| ArrowError::ExternalError(e.into()))?
                .map(|i| i.map_err(|e| ArrowError::ExternalError(e.into())));

            return write_ipc_err(query, batch_builder);
        }
        let records = self
            .reader
            .records(&self.header)
            .map(|i| i.map_err(|e| ArrowError::ExternalError(e.into())));
        write_ipc_err(records, batch_builder)
    }

    pub fn records_to_ipc_from_vpos(
        &mut self,
        pos_lo: (u64, u16),
        pos_hi: (u64, u16),
        tags: Option<HashSet<&str>>,
    ) -> Result<Vec<u8>, ArrowError> {
        let vpos_lo = bgzf::VirtualPosition::try_from(pos_lo)
            .map_err(|e| ArrowError::ExternalError(e.into()))?;
        let vpos_hi = bgzf::VirtualPosition::try_from(pos_hi)
            .map_err(|e| ArrowError::ExternalError(e.into()))?;
        let batch_builder = BamBatchBuilder::new(1024, &self.header, tags)?;
        let records = BamRecords::new(&mut self.reader, &self.header, vpos_lo, vpos_hi)
            .map(|i| i.map_err(|e| ArrowError::ExternalError(e.into())));
        write_ipc_err(records, batch_builder)
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
    tags: Option<HashSet<&'a str>>,
    tag_values: HashMap<String, GenericStringBuilder<i32>>,
}

impl<'a> BamBatchBuilder<'a> {
    pub fn new(
        capacity: usize,
        header: &'a sam::Header,
        tags: Option<HashSet<&'a str>>,
    ) -> Result<Self, ArrowError> {
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
            tags,
            tag_values: HashMap::new(),
        })
    }

    fn add_tags(&mut self, record: &sam::alignment::Record) {
        // Add tags from a record to the tag_values hashmap. This is a method
        // That is called by the BatchBuilder trait implementation but it's not
        // part of the trait itself.
        let tags;
        match &self.tags {
            Some(t) => {
                tags = t
                    .iter()
                    .filter_map(|x| Tag::from_str(x).ok())
                    .collect::<HashSet<_>>()
            }

            // if no tags are specified, return all tags
            None => tags = record.data().keys().collect::<HashSet<_>>(),
        }

        // Go through each expected tag (the ones asked for or all tags if the asked for were omitted)
        for tag in tags {
            let tag_str = tag.to_string();

            match record.data().get(&tag) {
                // Does the record have the tag we're looking for?
                // If it doesn't, we don't care and move on.
                Some(value) => match self.tag_values.get_mut(tag_str.as_str()) {
                    // See if we're already keeping track of this tag
                    Some(tag_value) => {
                        // If we are, we need to add in nulls for all the previous records
                        // that didn't have this tag so that in the end all returned columns
                        // have the same length
                        while tag_value.len() < self.qname.len() - 1 {
                            tag_value.append_null();
                        }

                        // Add the actual value
                        tag_value.append_value(value.to_string());
                    }
                    None => {
                        // We haven't seen this tag before so we have to create a GenericStringBuilder
                        // to start adding values.
                        let mut str_builder = GenericStringBuilder::<i32>::new();

                        // Add null values for all the previous rows where this tag wasn't present
                        while str_builder.len() < self.qname.len() - 1 {
                            str_builder.append_null();
                        }

                        // Add the actual value
                        str_builder.append_value(value.to_string());

                        // Attach this GenericStringBuilder to the tag
                        self.tag_values.insert(tag_str, str_builder);
                    }
                },
                None => {}
            }
        }
    }
}

impl<'a> BatchBuilder for BamBatchBuilder<'a> {
    type Record<'x> = &'x sam::alignment::Record;

    fn push(&mut self, record: Self::Record<'_>) {
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

        // Add tag values. This is a little more complicated than the other fields so we're
        // doing it in a separate method.
        self.add_tags(record);
    }

    fn finish(mut self) -> Result<RecordBatch, ArrowError> {
        let num_records = self.qname.len();

        let mut record_batch = vec![
            // spec
            (
                String::from("qname"),
                Arc::new(self.qname.finish()) as ArrayRef,
            ),
            (
                String::from("flag"),
                Arc::new(self.flag.finish()) as ArrayRef,
            ),
            (
                String::from("rname"),
                Arc::new(self.rname.finish()) as ArrayRef,
            ),
            (String::from("pos"), Arc::new(self.pos.finish()) as ArrayRef),
            (
                String::from("mapq"),
                Arc::new(self.mapq.finish()) as ArrayRef,
            ),
            (
                String::from("cigar"),
                Arc::new(self.cigar.finish()) as ArrayRef,
            ),
            (
                String::from("rnext"),
                Arc::new(self.rnext.finish()) as ArrayRef,
            ),
            (
                String::from("pnext"),
                Arc::new(self.pnext.finish()) as ArrayRef,
            ),
            (
                String::from("tlen"),
                Arc::new(self.tlen.finish()) as ArrayRef,
            ),
            (String::from("seq"), Arc::new(self.seq.finish()) as ArrayRef),
            (
                String::from("qual"),
                Arc::new(self.qual.finish()) as ArrayRef,
            ),
            // extra
            (String::from("end"), Arc::new(self.end.finish()) as ArrayRef),
        ];

        // Add the tags
        for (key, value) in self.tag_values.iter_mut() {
            // If a tag appeared in some records but not the last ones, we need to
            // fill in the remaining blanks with nulls
            while value.len() < num_records {
                value.append_null();
            }

            let array_ref = Arc::new(value.finish()) as ArrayRef;

            record_batch.push((key.to_string(), array_ref));
        }

        RecordBatch::try_from_iter(record_batch)
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
        vpos_hi: bgzf::VirtualPosition,
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

    fn read_record_batch(region: Option<&str>, tags: Option<HashSet<&str>>) -> RecordBatch {
        let mut dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push("../fixtures/sample.bam");
        let mut reader = BamReader::new_from_path(dir.to_str().unwrap()).unwrap();
        let ipc = reader.records_to_ipc(region, tags).unwrap();
        let cursor = std::io::Cursor::new(ipc);
        let mut arrow_reader = FileReader::try_new(cursor, None).unwrap();
        // make sure we have one batch
        assert_eq!(arrow_reader.num_batches(), 1);
        arrow_reader.next().unwrap().unwrap()
    }

    #[test]
    fn test_read_tags() {
        let record_batch = read_record_batch(None, Some(HashSet::from(["MD"])));

        // Check to make sure the tags we requested are present
        assert!(record_batch.column_by_name("MD").is_some());

        // Check to make sure tags we didn't request are not present
        assert!(record_batch.column_by_name("NM").is_none());

        let record_batch = read_record_batch(None, None);

        // Check to make sure both tags are present when no tags are specified
        dbg!(&record_batch);

        assert!(record_batch.column_by_name("MD").is_some());
        assert!(record_batch.column_by_name("NM").is_some());

        assert_eq!(record_batch.num_rows(), 6);
    }

    #[test]
    fn test_read_all() {
        let record_batch = read_record_batch(None, None);
        assert_eq!(record_batch.num_rows(), 6);
    }

    #[test]
    fn test_region_full() {
        let record_batch = read_record_batch(Some("chr1"), None);
        assert_eq!(record_batch.num_rows(), 4);
    }

    #[test]
    fn rest_region_partial() {
        let record_batch = read_record_batch(Some("chr1:1-100000"), None);
        assert_eq!(record_batch.num_rows(), 2);
    }
}
