use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufReader, Read, Seek};
use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, Float32Builder, GenericStringBuilder, Int16Builder, Int32Array, Int32Builder,
    Int8Builder, StringArray, StringDictionaryBuilder, StructArray, UInt16Array, UInt16Builder,
    UInt32Builder, UInt8Array, UInt8Builder,
};
use arrow::{datatypes::Int32Type, error::ArrowError, record_batch::RecordBatch};
use noodles::core::Region;
use noodles::sam::record::data::field::Tag;
use noodles::sam::record::Data;
use noodles::{bam, bgzf, csi, sam};

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
    /// If the region is `None`, all records are returned.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use oxbow::bam::BamReader;
    ///
    /// let mut reader = BamReader::new_from_path("sample.bam").unwrap();
    /// let ipc = reader.records_to_ipc(Some("sq0:1-1000")).unwrap();
    /// ```
    pub fn records_to_ipc(&mut self, region: Option<&str>) -> Result<Vec<u8>, ArrowError> {
        let batch_builder = BamBatchBuilder::new(1024, &self.header)?;
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
    ) -> Result<Vec<u8>, ArrowError> {
        let vpos_lo = bgzf::VirtualPosition::try_from(pos_lo)
            .map_err(|e| ArrowError::ExternalError(e.into()))?;
        let vpos_hi = bgzf::VirtualPosition::try_from(pos_hi)
            .map_err(|e| ArrowError::ExternalError(e.into()))?;
        let batch_builder = BamBatchBuilder::new(1024, &self.header)?;
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
    tags: TagsBuilder,
}

enum TagArrayBuilder {
    Character(GenericStringBuilder<i32>),
    Int8(Int8Builder),
    UInt8(UInt8Builder),
    Int16(Int16Builder),
    UInt16(UInt16Builder),
    Int32(Int32Builder),
    UInt32(UInt32Builder),
    Float(Float32Builder),
    String(GenericStringBuilder<i32>),
    Hex(GenericStringBuilder<i32>),
}

struct TagsBuilder {
    inner: HashMap<Tag, TagArrayBuilder>,
    seen: usize,
    tags: HashSet<Tag>,
}

impl TagsBuilder {
    pub fn new() -> Self {
        Self {
            inner: HashMap::new(),
            seen: 0,
            tags: HashSet::new(),
        }
    }

    pub fn push_tags(&mut self, data: &'_ Data) {
        use sam::record::data::field::Value;
        self.tags.extend(data.keys());

        for tag in self.tags.iter() {
            match data.get(tag) {
                Some(Value::Character(v)) => {
                    let builder = self.inner.entry(*tag).or_insert_with(|| {
                        let mut builder = GenericStringBuilder::<i32>::new();
                        builder.extend(std::iter::repeat(None::<&str>).take(self.seen));
                        TagArrayBuilder::Character(builder)
                    });
                    match builder {
                        TagArrayBuilder::Character(builder) => {
                            builder.append_value(v.to_string());
                        }
                        _ => panic!("Wrong type"),
                    }
                }
                Some(Value::Int8(v)) => {
                    let builder = self.inner.entry(*tag).or_insert_with(|| {
                        let mut builder = Int8Builder::new();
                        builder.extend(std::iter::repeat(None).take(self.seen));
                        TagArrayBuilder::Int8(builder)
                    });
                    match builder {
                        TagArrayBuilder::Int8(builder) => {
                            builder.append_value(*v);
                        }
                        _ => panic!("Wrong type"),
                    }
                }
                Some(Value::UInt8(v)) => {
                    let builder = self.inner.entry(*tag).or_insert_with(|| {
                        let mut builder = UInt8Builder::new();
                        builder.extend(std::iter::repeat(None).take(self.seen));
                        TagArrayBuilder::UInt8(builder)
                    });
                    match builder {
                        TagArrayBuilder::UInt8(builder) => {
                            builder.append_value(*v);
                        }
                        _ => panic!("Wrong type"),
                    }
                }
                Some(Value::Int16(v)) => {
                    let builder = self.inner.entry(*tag).or_insert_with(|| {
                        let mut builder = Int16Builder::new();
                        builder.extend(std::iter::repeat(None).take(self.seen));
                        TagArrayBuilder::Int16(builder)
                    });
                    match builder {
                        TagArrayBuilder::Int16(builder) => {
                            builder.append_value(*v);
                        }
                        _ => panic!("Wrong type"),
                    }
                }
                Some(Value::UInt16(v)) => {
                    let builder = self.inner.entry(*tag).or_insert_with(|| {
                        let mut builder = UInt16Builder::new();
                        builder.extend(std::iter::repeat(None).take(self.seen));
                        TagArrayBuilder::UInt16(builder)
                    });
                    match builder {
                        TagArrayBuilder::UInt16(builder) => {
                            builder.append_value(*v);
                        }
                        _ => panic!("Wrong type"),
                    }
                }
                Some(Value::Int32(v)) => {
                    let builder = self.inner.entry(*tag).or_insert_with(|| {
                        let mut builder = Int32Builder::new();
                        builder.extend(std::iter::repeat(None).take(self.seen));
                        TagArrayBuilder::Int32(builder)
                    });
                    match builder {
                        TagArrayBuilder::Int32(builder) => {
                            builder.append_value(*v);
                        }
                        _ => panic!("Wrong type"),
                    }
                }
                Some(Value::UInt32(v)) => {
                    let builder = self.inner.entry(*tag).or_insert_with(|| {
                        let mut builder = UInt32Builder::new();
                        builder.extend(std::iter::repeat(None).take(self.seen));
                        TagArrayBuilder::UInt32(builder)
                    });
                    match builder {
                        TagArrayBuilder::UInt32(builder) => {
                            builder.append_value(*v);
                        }
                        _ => panic!("Wrong type"),
                    }
                }
                Some(Value::Float(v)) => {
                    let builder = self.inner.entry(*tag).or_insert_with(|| {
                        let mut builder = Float32Builder::new();
                        builder.extend(std::iter::repeat(None).take(self.seen));
                        TagArrayBuilder::Float(builder)
                    });
                    match builder {
                        TagArrayBuilder::Float(builder) => {
                            builder.append_value(*v);
                        }
                        _ => panic!("Wrong type"),
                    }
                }
                Some(Value::String(v)) => {
                    let builder = self.inner.entry(*tag).or_insert_with(|| {
                        let mut builder = GenericStringBuilder::<i32>::new();
                        builder.extend(std::iter::repeat(None::<&str>).take(self.seen));
                        TagArrayBuilder::String(builder)
                    });
                    match builder {
                        TagArrayBuilder::String(builder) => {
                            builder.append_value(v.as_str());
                        }
                        _ => panic!("Wrong type"),
                    }
                }
                Some(Value::Hex(v)) => {
                    let builder = self.inner.entry(*tag).or_insert_with(|| {
                        let mut builder = GenericStringBuilder::<i32>::new();
                        builder.extend(std::iter::repeat(None::<&str>).take(self.seen));
                        TagArrayBuilder::Hex(builder)
                    });
                    match builder {
                        TagArrayBuilder::Hex(builder) => {
                            builder.append_value(v.as_ref());
                        }
                        _ => panic!("Wrong type"),
                    }
                }
                Some(Value::Array(_)) => {
                    panic!("Array not implemented");
                }
                None => match self.inner.get_mut(tag).expect("Tag not found") {
                    TagArrayBuilder::Character(builder) => builder.append_null(),
                    TagArrayBuilder::Int8(builder) => builder.append_null(),
                    TagArrayBuilder::UInt8(builder) => builder.append_null(),
                    TagArrayBuilder::Int16(builder) => builder.append_null(),
                    TagArrayBuilder::UInt16(builder) => builder.append_null(),
                    TagArrayBuilder::Int32(builder) => builder.append_null(),
                    TagArrayBuilder::UInt32(builder) => builder.append_null(),
                    TagArrayBuilder::Float(builder) => builder.append_null(),
                    TagArrayBuilder::String(builder) => builder.append_null(),
                    TagArrayBuilder::Hex(builder) => builder.append_null(),
                },
            }
        }
        self.seen += 1;
    }

    fn try_finish(&mut self) -> Result<StructArray, ArrowError> {
        let keys = self.inner.keys().map(|x| x.to_string()).collect::<Vec<_>>();
        let arrays: Vec<(&str, ArrayRef)> = self
            .inner
            .iter_mut()
            .zip(keys.iter())
            .map(|((_, builder), key)| {
                let arr: ArrayRef = match builder {
                    TagArrayBuilder::Character(builder) => Arc::new(builder.finish()) as ArrayRef,
                    TagArrayBuilder::Int8(builder) => Arc::new(builder.finish()) as ArrayRef,
                    TagArrayBuilder::UInt8(builder) => Arc::new(builder.finish()) as ArrayRef,
                    TagArrayBuilder::Int16(builder) => Arc::new(builder.finish()) as ArrayRef,
                    TagArrayBuilder::UInt16(builder) => Arc::new(builder.finish()) as ArrayRef,
                    TagArrayBuilder::Int32(builder) => Arc::new(builder.finish()) as ArrayRef,
                    TagArrayBuilder::UInt32(builder) => Arc::new(builder.finish()) as ArrayRef,
                    TagArrayBuilder::Float(builder) => Arc::new(builder.finish()) as ArrayRef,
                    TagArrayBuilder::String(builder) => Arc::new(builder.finish()) as ArrayRef,
                    TagArrayBuilder::Hex(builder) => Arc::new(builder.finish()) as ArrayRef,
                };
                (key.as_str(), arr)
            })
            .collect();

        StructArray::try_from(arrays)
    }
}

impl Default for TagsBuilder {
    fn default() -> Self {
        Self::new()
    }
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
            tags: TagsBuilder::new(),
        })
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
        self.tags.push_tags(record.data());
    }

    fn finish(mut self) -> Result<RecordBatch, ArrowError> {
        let tags = self.tags.try_finish()?;
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
            ("tags", Arc::new(tags) as ArrayRef),
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

    fn read_record_batch(region: Option<&str>) -> RecordBatch {
        let mut dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push("../fixtures/sample.bam");
        let mut reader = BamReader::new_from_path(dir.to_str().unwrap()).unwrap();
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
