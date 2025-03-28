use std::io;
use std::sync::Arc;

use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use arrow::record_batch::RecordBatchReader;

use noodles::core::region::Interval;
use noodles::csi::binning_index;

use crate::bed::model::BatchBuilder;
use crate::bed::model::Push as _;
use crate::util::query::BgzfChunkReader;

/// A record batch iterator yielding BED records that intersect a genomic range.
pub struct BatchIterator<R> {
    reader: R,
    builder: BatchBuilder,
    batch_size: usize,
    limit: usize,
    count: usize,
    header: binning_index::index::Header,
    reference_sequence_id: usize,
    interval: Interval,
}

impl<R> BatchIterator<R> {
    pub fn new(
        reader: R,
        header: binning_index::index::Header,
        reference_sequence_id: usize,
        interval: Interval,
        builder: BatchBuilder,
        batch_size: usize,
        limit: Option<usize>,
    ) -> Self {
        Self {
            reader,
            builder,
            batch_size,
            limit: limit.unwrap_or(usize::MAX),
            count: 0,
            reference_sequence_id,
            interval,
            header,
        }
    }
}

impl<R> RecordBatchReader for BatchIterator<R>
where
    BatchIterator<R>: Iterator<Item = Result<RecordBatch, ArrowError>>,
{
    fn schema(&self) -> arrow::datatypes::SchemaRef {
        Arc::new(self.builder.get_arrow_schema())
    }
}

impl<R> Iterator for BatchIterator<noodles::bed::io::Reader<3, BgzfChunkReader<R>>>
where
    R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = noodles::bed::Record::default();
        let mut count = 0;

        while count < self.batch_size && self.count < self.limit {
            match self.reader.read_record(&mut record) {
                Ok(0) => break,
                Ok(_) => {
                    match intersects(
                        &self.header,
                        &record,
                        self.reference_sequence_id,
                        self.interval,
                    ) {
                        Ok(true) => match self.builder.push(&record) {
                            Ok(()) => {
                                self.count += 1;
                                count += 1;
                            }
                            Err(e) => return Some(Err(e.into())),
                        },
                        Ok(false) => {}
                        Err(e) => return Some(Err(e.into())),
                    }
                }
                Err(e) => return Some(Err(e.into())),
            };
        }

        if count == 0 {
            None
        } else {
            let batch = self.builder.finish();
            Some(batch)
        }
    }
}

pub fn intersects(
    header: &noodles::csi::binning_index::index::Header,
    record: &noodles::bed::Record<3>,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> io::Result<bool> {
    let rname = record.reference_sequence_name();
    match (
        header.reference_sequence_names().get_index_of(rname),
        record.feature_start(),
        record.feature_end(),
    ) {
        (Some(id), Ok(start), Some(Ok(end))) => {
            let alignment_interval = (start..=end).into();
            Ok(id == reference_sequence_id && region_interval.intersects(alignment_interval))
        }
        _ => Ok(false),
    }
}
