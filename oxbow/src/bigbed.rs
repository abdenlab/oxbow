use arrow::array::StringArray;
use arrow::array::{
    ArrayRef, Float32Array, Float32Builder, Float64Array, Float64Builder, GenericByteBuilder,
    Int16Array, Int16Builder, Int32Array, Int32Builder, Int8Array, Int8Builder, StringBuilder,
    StringDictionaryBuilder, UInt16Array, UInt16Builder, UInt32Array, UInt32Builder, UInt8Array,
    UInt8Builder,
};
use arrow::datatypes::Int32Type;
use arrow::{error::ArrowError, record_batch::RecordBatch};
use bigtools::utils::reopen::ReopenableFile;
use bigtools::{BBIRead, BigBedRead};
use noodles::core::Region;
use std::collections::HashSet;
use std::io::{Read, Seek};
use std::sync::Arc;

use crate::batch_builder::{finish_batch, BatchBuilder};

/// A BigBed reader.
pub struct BigBedReader<R> {
    read: BigBedRead<R>,
}

pub struct BigBedRecord<'a> {
    pub chrom: &'a str,
    pub start: u32,
    pub end: u32,
    pub rest: &'a str,
}

impl BigBedReader<ReopenableFile> {
    pub fn new_from_path(path: &str) -> std::io::Result<Self> {
        let read = BigBedRead::open_file(path)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
        Ok(Self { read })
    }
}
impl<R: Read + Seek> BigBedReader<R> {
    /// Creates a BigBed reader.
    pub fn new(read: R) -> std::io::Result<Self> {
        let read = BigBedRead::open(read)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
        Ok(Self { read })
    }

    /// Returns the records in the given region as Apache Arrow IPC.
    ///
    /// If the region is `None`, all records are returned.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use oxbow::bigbed::BigBedReader;
    ///
    /// let mut reader = BigBedReader::new_from_path("sample.bigBed").unwrap();
    /// let ipc = reader.records_to_ipc(Some("sq0:1-1000"), None).unwrap();
    /// ```
    pub fn records_to_ipc(
        &mut self,
        region: Option<&str>,
        fields: Option<HashSet<&str>>,
    ) -> Result<Vec<u8>, ArrowError> {
        let mut batch_builder = BigBedBatchBuilder::new(1024, &mut self.read, fields)?;
        match region {
            Some(region) => {
                let region: Region = region.parse().unwrap();
                let chrom_name = region.name().to_owned();
                let (start, end) = match (region.interval().start(), region.interval().end()) {
                    (Some(start), Some(end)) => {
                        let start = start.get() as u32 - 1; // 1-based to 0-based
                        let end = end.get() as u32;
                        (start, end)
                    }
                    (Some(start), None) => {
                        let start = start.get() as u32 - 1; // 1-based to 0-based
                        let end = self
                            .read
                            .get_chroms()
                            .iter()
                            .find(|c| c.name == chrom_name)
                            .map(|c| c.length);
                        let end = end.ok_or_else(|| {
                            ArrowError::InvalidArgumentError("Invalid chromosome".to_string())
                        })?;
                        (start, end)
                    }
                    (None, Some(end)) => {
                        let start = 0;
                        let end = end.get() as u32;
                        (start, end)
                    }
                    (None, None) => {
                        let start = 0;
                        let end = self
                            .read
                            .get_chroms()
                            .iter()
                            .find(|c| c.name == chrom_name)
                            .map(|c| c.length);
                        let end = end.ok_or_else(|| {
                            ArrowError::InvalidArgumentError("Invalid chromosome".to_string())
                        })?;
                        (start, end)
                    }
                };
                let values = match self.read.get_interval(&chrom_name, start, end) {
                    Ok(v) => v,
                    Err(e) => {
                        return Err(ArrowError::ExternalError(Box::new(e)));
                    }
                };
                for value in values {
                    let v = value.unwrap();
                    let record = BigBedRecord {
                        chrom: &chrom_name,
                        start: v.start,
                        end: v.end,
                        rest: &v.rest,
                    };
                    batch_builder.push(record);
                }
                finish_batch(batch_builder)
            }
            None => {
                // Can't use write_ipc, because we have separate iterators for each chrom
                let chroms = self.read.get_chroms().into_iter();
                for chrom in chroms {
                    let start = 0;
                    let end = chrom.length;
                    let values = self.read.get_interval(&chrom.name, start, end).unwrap();
                    for value in values {
                        let v = value.unwrap();
                        let record = BigBedRecord {
                            chrom: &chrom.name,
                            start: v.start,
                            end: v.end,
                            rest: &v.rest,
                        };
                        batch_builder.push(record);
                    }
                }
                finish_batch(batch_builder)
            }
        }
    }
}

enum Column {
    Int(Int32Builder),
    Uint(UInt32Builder),
    Short(Int16Builder),
    Ushort(UInt16Builder),
    Byte(Int8Builder),
    Ubyte(UInt8Builder),
    Float(Float32Builder),
    Double(Float64Builder),
    String(StringBuilder),
}

enum ExtraColumns {
    Full(StringBuilder),
    Split(Vec<Option<(String, Column)>>),
}

struct BigBedBatchBuilder {
    chrom: StringDictionaryBuilder<Int32Type>,
    start: UInt32Builder,
    end: UInt32Builder,
    extra: ExtraColumns,
}

impl BigBedBatchBuilder {
    pub fn new<R: Read + Seek>(
        capacity: usize,
        read: &mut BigBedRead<R>,
        columns: Option<HashSet<&str>>,
    ) -> Result<Self, ArrowError> {
        let chroms: Vec<String> = read.get_chroms().iter().map(|c| c.name.clone()).collect();
        let chroms: StringArray = StringArray::from(chroms);
        let autosql = read.autosql().unwrap();
        let mut declarations = bigtools::bed::autosql::parse::parse_autosql(&autosql).unwrap();
        if declarations.len() > 1 {
            panic!("Unexpected extra declarations");
        }
        let declaration = declarations.pop();
        let extra = match declaration {
            None => ExtraColumns::Full(GenericByteBuilder::<_>::with_capacity(capacity, 0)),
            Some(declaration) => {
                // First 3 are chrom, start, end
                let columns = declaration
                    .fields
                    .into_iter()
                    .skip(3)
                    .map(|f| {
                        if columns
                            .as_ref()
                            .map_or(false, |c| !c.contains(f.name.as_str()))
                        {
                            return None;
                        }
                        use bigtools::bed::autosql::parse::FieldType::*;
                        let column = match f.field_type {
                            Int => Column::Int(Int32Array::builder(capacity)),
                            Uint => Column::Uint(UInt32Array::builder(capacity)),
                            Short => Column::Short(Int16Array::builder(capacity)),
                            Ushort => Column::Ushort(UInt16Array::builder(capacity)),
                            Byte => Column::Byte(Int8Array::builder(capacity)),
                            Ubyte => Column::Ubyte(UInt8Array::builder(capacity)),
                            Float => Column::Float(Float32Array::builder(capacity)),
                            Double => Column::Double(Float64Array::builder(capacity)),
                            Char => {
                                Column::String(GenericByteBuilder::<_>::with_capacity(capacity, 1))
                            }
                            String => {
                                Column::String(GenericByteBuilder::<_>::with_capacity(capacity, 0))
                            }
                            Lstring => {
                                Column::String(GenericByteBuilder::<_>::with_capacity(capacity, 0))
                            }
                            Bigint => Column::Int(Int32Array::builder(capacity)),
                            Enum(_) => {
                                Column::String(GenericByteBuilder::<_>::with_capacity(capacity, 1))
                            }
                            Set(_) => {
                                Column::String(GenericByteBuilder::<_>::with_capacity(capacity, 1))
                            }
                            Declaration(_, _) => panic!("Unexpected field type"),
                        };
                        Some((f.name, column))
                    })
                    .collect();
                ExtraColumns::Split(columns)
            }
        };
        Ok(Self {
            chrom: StringDictionaryBuilder::<Int32Type>::new_with_dictionary(capacity, &chroms)?,
            start: UInt32Array::builder(capacity),
            end: UInt32Array::builder(capacity),
            extra,
        })
    }
}

impl BatchBuilder for BigBedBatchBuilder {
    type Record<'a> = BigBedRecord<'a>;

    fn push(&mut self, record: Self::Record<'_>) {
        self.chrom.append_value(record.chrom);
        self.start.append_value(record.start);
        self.end.append_value(record.end);
        match &mut self.extra {
            ExtraColumns::Full(rest) => rest.append_value(record.rest),
            ExtraColumns::Split(columns) => {
                for (builder, col) in
                    std::iter::zip(columns.iter_mut(), record.rest.split_whitespace())
                {
                    let Some((_, builder)) = builder else {
                        continue;
                    };
                    match builder {
                        Column::Int(builder) => {
                            let value: i32 = col.replace(",", "").parse().unwrap();
                            builder.append_value(value);
                        }
                        Column::Uint(builder) => {
                            let value: u32 = col.replace(",", "").parse().unwrap();
                            builder.append_value(value);
                        }
                        Column::Short(builder) => {
                            let value: i16 = col.replace(",", "").parse().unwrap();
                            builder.append_value(value);
                        }
                        Column::Ushort(builder) => {
                            let value: u16 = col.replace(",", "").parse().unwrap();
                            builder.append_value(value);
                        }
                        Column::Byte(builder) => {
                            let value: i8 = col.replace(",", "").parse().unwrap();
                            builder.append_value(value);
                        }
                        Column::Ubyte(builder) => {
                            let value: u8 = col.replace(",", "").parse().unwrap();
                            builder.append_value(value);
                        }
                        Column::Float(builder) => {
                            let value: f32 = col.replace(",", "").parse().unwrap();
                            builder.append_value(value);
                        }
                        Column::Double(builder) => {
                            let value: f64 = col.replace(",", "").parse().unwrap();
                            builder.append_value(value);
                        }
                        Column::String(builder) => {
                            builder.append_value(col);
                        }
                    }
                }
            }
        }
    }

    fn finish(mut self) -> Result<RecordBatch, ArrowError> {
        match self.extra {
            ExtraColumns::Full(mut builder) => RecordBatch::try_from_iter(vec![
                ("chrom", Arc::new(self.chrom.finish()) as ArrayRef),
                ("start", Arc::new(self.start.finish()) as ArrayRef),
                ("end", Arc::new(self.end.finish()) as ArrayRef),
                ("rest", Arc::new(builder.finish()) as ArrayRef),
            ]),
            ExtraColumns::Split(columns) => RecordBatch::try_from_iter(
                vec![
                    (
                        "chrom".to_string(),
                        Arc::new(self.chrom.finish()) as ArrayRef,
                    ),
                    (
                        "start".to_string(),
                        Arc::new(self.start.finish()) as ArrayRef,
                    ),
                    ("end".to_string(), Arc::new(self.end.finish()) as ArrayRef),
                ]
                .into_iter()
                .chain(columns.into_iter().flatten().map(|mut c| {
                    let array = match &mut c.1 {
                        Column::Int(builder) => Arc::new(builder.finish()) as ArrayRef,
                        Column::Uint(builder) => Arc::new(builder.finish()) as ArrayRef,
                        Column::Short(builder) => Arc::new(builder.finish()) as ArrayRef,
                        Column::Ushort(builder) => Arc::new(builder.finish()) as ArrayRef,
                        Column::Byte(builder) => Arc::new(builder.finish()) as ArrayRef,
                        Column::Ubyte(builder) => Arc::new(builder.finish()) as ArrayRef,
                        Column::Float(builder) => Arc::new(builder.finish()) as ArrayRef,
                        Column::Double(builder) => Arc::new(builder.finish()) as ArrayRef,
                        Column::String(builder) => Arc::new(builder.finish()) as ArrayRef,
                    };
                    (c.0, array)
                })),
            ),
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
        dir.push("../fixtures/small.bigBed");
        let mut reader = BigBedReader::new_from_path(dir.to_str().unwrap()).unwrap();
        let ipc = reader.records_to_ipc(region, None).unwrap();
        let cursor = std::io::Cursor::new(ipc);
        let mut arrow_reader = FileReader::try_new(cursor, None).unwrap();
        // make sure we have one batch
        assert_eq!(arrow_reader.num_batches(), 1);
        arrow_reader.next().unwrap().unwrap()
    }
    /*
        #[test]
        fn test_read_all() {
            let record_batch = read_record_batch(None);
            assert_eq!(record_batch.num_rows(), 100000);
        }

        #[test]
        fn test_region_full() {
            let record_batch = read_record_batch(Some("chr17"));
            assert_eq!(record_batch.num_rows(), 37385);
        }
    */
    #[test]
    fn rest_region_partial() {
        let record_batch = read_record_batch(Some("chr17:100000-150000"));
        assert_eq!(record_batch.num_rows(), 3);
    }
}
