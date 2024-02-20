use arrow::{
    array::{
        ArrayRef, GenericStringBuilder, ListBuilder, StructBuilder, UInt16Builder, UInt64Builder,
        UInt8Builder,
    },
    datatypes::{DataType, Field},
    error::ArrowError,
    record_batch::RecordBatch,
};
use noodles::bed;
use std::{
    fs::File,
    io::{self, BufRead, BufReader},
    str,
    sync::Arc,
};

use crate::batch_builder::{write_ipc, BatchBuilder};

pub fn new_from_path(path: &str) -> io::Result<bed::Reader<BufReader<File>>> {
    File::open(path).map(BufReader::new).map(bed::Reader::new)
}

pub fn new_from_reader<R>(bed: R) -> io::Result<bed::Reader<BufReader<R>>>
where
    R: io::Read,
{
    Ok(bed::Reader::new(BufReader::new(bed)))
}

/// Returns the records in the given region as Apache Arrow IPC.
///
/// If the region is `None`, all records are returned.
///
/// # Examples
///
/// ```no_run
/// use oxbow::bed;
///
/// let mut reader = bed::new_from_path("sample.bed").unwrap();
/// let ipc = bed::records_to_ipc(reader).unwrap();
/// ```
pub fn records_to_ipc<T>(mut reader: bed::Reader<T>) -> Result<Vec<u8>, ArrowError>
where
    T: BufRead,
{
    let batch_builder = BedBatchBuilder::new(1024)?;
    let records = reader.records().map(|r| r.unwrap());
    write_ipc(records, batch_builder)
}

struct BedBatchBuilder {
    reference_sequence_name: GenericStringBuilder<i32>,
    start_position: UInt64Builder,
    end_position: UInt64Builder,
    name: GenericStringBuilder<i32>,
    score: UInt16Builder,
    strand: GenericStringBuilder<i32>,
    thick_start: UInt64Builder,
    thick_end: UInt64Builder,
    color: StructBuilder,
    block_sizes: ListBuilder<UInt64Builder>,
    block_starts: ListBuilder<UInt64Builder>,
}

impl BedBatchBuilder {
    pub fn new(_capacity: usize) -> Result<Self, ArrowError> {
        Ok(Self {
            reference_sequence_name: GenericStringBuilder::<i32>::new(),
            start_position: UInt64Builder::new(),
            end_position: UInt64Builder::new(),
            name: GenericStringBuilder::<i32>::new(),
            score: UInt16Builder::new(),
            strand: GenericStringBuilder::<i32>::new(),
            thick_start: UInt64Builder::new(),
            thick_end: UInt64Builder::new(),
            color: StructBuilder::new(
                vec![
                    Field::new("r", DataType::UInt8, true),
                    Field::new("g", DataType::UInt8, true),
                    Field::new("b", DataType::UInt8, true),
                ],
                vec![
                    Box::new(UInt8Builder::new()),
                    Box::new(UInt8Builder::new()),
                    Box::new(UInt8Builder::new()),
                ],
            ),
            block_sizes: ListBuilder::new(UInt64Builder::new()),
            block_starts: ListBuilder::new(UInt64Builder::new()),
        })
    }
}

impl BatchBuilder for BedBatchBuilder {
    // Noodles implements bed Record structs for bed3-bed9 and bed12.
    // The following only handles bed12.
    type Record<'a> = &'a bed::Record<12>;

    fn push(&mut self, record: Self::Record<'_>) {
        // Mandatory fields (spec)
        self.reference_sequence_name
            .append_value(record.reference_sequence_name());
        self.start_position
            .append_value(record.start_position().get().try_into().unwrap());
        self.end_position
            .append_value(record.end_position().get().try_into().unwrap());

        // Optional fields (spec)
        if let Some(name) = record.name() {
            self.name.append_value(name.to_string());
        } else {
            self.name.append_null();
        }

        if let Some(score) = record.score() {
            self.score.append_value(score.get().try_into().unwrap());
        } else {
            // A score of 0 will be handled as a null as well
            self.score.append_null();
        }

        if let Some(strand) = record.strand() {
            self.strand.append_value(strand.to_string());
        } else {
            self.strand.append_null();
        }

        // thick_start and thick_end are required fields in Noodle's bed7/8 structs.
        self.thick_start
            .append_value(record.thick_start().get().try_into().unwrap());

        self.thick_end
            .append_value(record.thick_end().get().try_into().unwrap());

        if let Some(color) = record.color() {
            let colors = [color.red(), color.green(), color.blue()];
            for (i, val) in colors.iter().enumerate() {
                self.color
                    .field_builder::<UInt8Builder>(i)
                    .unwrap()
                    .append_value(*val);
            }

            // Handle 1 or 2 missing colors
            if colors.len() == 2 {
                self.color
                    .field_builder::<UInt8Builder>(2)
                    .unwrap()
                    .append_null();
            } else if colors.len() == 1 {
                self.color
                    .field_builder::<UInt8Builder>(1)
                    .unwrap()
                    .append_null();
                self.color
                    .field_builder::<UInt8Builder>(2)
                    .unwrap()
                    .append_null();
            }
            self.color.append(true);
        } else {
            for i in 0..3 {
                self.color
                    .field_builder::<UInt8Builder>(i)
                    .unwrap()
                    .append_null();
            }
            self.color.append(true);
        }

        let (block_sizes, block_starts): (Vec<_>, Vec<_>) = record.blocks().iter().cloned().unzip();
        let block_sizes_option_u64: Vec<Option<u64>> = block_sizes
            .into_iter()
            .map(|usize_value| Some(usize_value as u64))
            .collect();
        let block_starts_option_u64: Vec<Option<u64>> = block_starts
            .into_iter()
            .map(|usize_value| Some(usize_value as u64))
            .collect();
        self.block_starts.append_value(block_sizes_option_u64);
        self.block_sizes.append_value(block_starts_option_u64);
    }

    fn finish(mut self) -> Result<RecordBatch, ArrowError> {
        RecordBatch::try_from_iter(vec![
            (
                "reference_sequence_name",
                Arc::new(self.reference_sequence_name.finish()) as ArrayRef,
            ),
            (
                "start_position",
                Arc::new(self.start_position.finish()) as ArrayRef,
            ),
            (
                "end_position",
                Arc::new(self.end_position.finish()) as ArrayRef,
            ),
            ("name", Arc::new(self.name.finish()) as ArrayRef),
            ("score", Arc::new(self.score.finish()) as ArrayRef),
            ("strand", Arc::new(self.strand.finish()) as ArrayRef),
            (
                "thick_start",
                Arc::new(self.thick_start.finish()) as ArrayRef,
            ),
            ("thick_end", Arc::new(self.thick_end.finish()) as ArrayRef),
            ("color", Arc::new(self.color.finish()) as ArrayRef),
            (
                "block_sizes",
                Arc::new(self.block_sizes.finish()) as ArrayRef,
            ),
            (
                "block_starts",
                Arc::new(self.block_starts.finish()) as ArrayRef,
            ),
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::ipc::reader::FileReader;
    use arrow::record_batch::RecordBatch;

    fn read_record_batch(fixture_path: &str) -> RecordBatch {
        let mut dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push(fixture_path);
        let reader = new_from_path(dir.to_str().unwrap()).unwrap();
        let ipc = records_to_ipc(reader).unwrap();
        let cursor = std::io::Cursor::new(ipc);
        let mut arrow_reader = FileReader::try_new(cursor, None).unwrap();
        // make sure we have one batch
        assert_eq!(arrow_reader.num_batches(), 1);
        arrow_reader.next().unwrap().unwrap()
    }

    //#[test]
    //fn test_read_all_bed3() {
    //    let record_batch = read_record_batch("../fixtures/bed3.bed");
    //    assert_eq!(record_batch.num_rows(), 2);
    //}

    //#[test]
    //fn test_read_all_bed4() {
    //    let record_batch = read_record_batch("../fixtures/bed4.bed");
    //    assert_eq!(record_batch.num_rows(), 2);
    //}

    //#[test]
    //fn test_read_all_bed5() {
    //    let record_batch = read_record_batch("../fixtures/bed5.bed");
    //    assert_eq!(record_batch.num_rows(), 2);
    //}

    //#[test]
    //fn test_read_all_bed6() {
    //    let record_batch = read_record_batch("../fixtures/bed6.bed");
    //    assert_eq!(record_batch.num_rows(), 2);
    //}

    //#[test]
    //fn test_read_all_bed7() {
    //    let record_batch = read_record_batch("../fixtures/bed7.bed");
    //    assert_eq!(record_batch.num_rows(), 2);
    //}

    //#[test]
    //fn test_read_all_bed8() {
    //    let record_batch = read_record_batch("../fixtures/bed8.bed");
    //    assert_eq!(record_batch.num_rows(), 2);
    //}

    //#[test]
    //fn test_read_all_bed9() {
    //    let record_batch = read_record_batch("../fixtures/bed9.bed");
    //    assert_eq!(record_batch.num_rows(), 2);
    //}

    #[test]
    fn test_read_all_bed12() {
        let record_batch = read_record_batch("../fixtures/bed12.bed");
        assert_eq!(record_batch.num_rows(), 2);
    }

    //#[test]
    //fn test_read_all_bed12plus() {
    //    let record_batch = read_record_batch("../fixtures/bed12plus.bed");
    //    assert_eq!(record_batch.num_rows(), 2);
    //}
}
