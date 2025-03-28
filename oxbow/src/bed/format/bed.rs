use std::io::{self, BufRead, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema;
use noodles::csi::BinningIndex;

use crate::bed::batch_iterator::{BatchIterator, QueryBatchIterator};
use crate::bed::model::BatchBuilder;
use crate::bed::model::BedSchema;
use crate::util::query::BgzfChunkReader;

/// A BED scanner.
///
/// # Examples
///
/// ```no_run
/// use oxbow::bed::format::bed::Scanner;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let inner = File::open("sample.bed").map(BufReader::new).unwrap();
/// let mut fmt_reader = noodles::bed::io::Reader::new(inner);
///
/// let bed_schema = "bed6+3".parse().unwrap();
/// let scanner = Scanner::new(bed_schema);
/// let batches = scanner.scan(fmt_reader, None, None, Some(1000)).unwrap();
/// ```
pub struct Scanner {
    bed_schema: BedSchema,
}

impl Scanner {
    /// Creates a BED scanner from a BED schema specifier.
    pub fn new(bed_schema: BedSchema) -> Self {
        Self { bed_schema }
    }

    /// Returns the standard field names.
    pub fn field_names(&self) -> Vec<String> {
        self.bed_schema.field_names()
    }

    /// Returns the Arrow schema.
    pub fn schema(&self, fields: Option<Vec<String>>) -> io::Result<Schema> {
        let batch_builder = BatchBuilder::new(fields, &self.bed_schema, 0)?;
        Ok(batch_builder.get_arrow_schema())
    }
}

impl Scanner {
    /// Returns an iterator over record batches.
    ///
    /// The scan will begin at the current position of the reader and will
    /// move the cursor to the end of the last record scanned.
    pub fn scan<R: BufRead>(
        &self,
        fmt_reader: noodles::bed::io::Reader<3, R>,
        fields: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = BatchBuilder::new(fields, &self.bed_schema, batch_size)?;
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }

    /// Returns an iterator over record batches satisfying a genomic range query.
    ///
    /// This operation requires a BGZF source and an Index.
    ///
    /// The scan will consume contiguous "chunks" of BGZF blocks and filter for
    /// records that overlap the given region. The cursor will move to the end
    /// of the last record scanned.
    pub fn scan_query<R: BufRead + Seek>(
        &self,
        fmt_reader: noodles::bed::io::Reader<3, noodles::bgzf::Reader<R>>,
        region: noodles::core::Region,
        index: impl BinningIndex,
        fields: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let reference_sequence_name = region.name().to_string();
        let interval = region.interval();

        let batch_builder = BatchBuilder::new(fields, &self.bed_schema, batch_size)?;

        let Some(header) = index.header() else {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                "Index header not found.",
            ));
        };
        let reference_sequence_id = resolve_chrom_id(header, &reference_sequence_name)?;
        let chunks = index.query(reference_sequence_id, interval)?;
        let bgzf_reader = fmt_reader.into_inner();
        let query_reader = BgzfChunkReader::new(bgzf_reader, chunks);
        let fmt_reader = noodles::bed::io::Reader::new(query_reader);
        let batch_iter = QueryBatchIterator::new(
            fmt_reader,
            header.clone(),
            reference_sequence_id,
            interval,
            batch_builder,
            batch_size,
            limit,
        );
        Ok(batch_iter)
    }
}

fn resolve_chrom_id(
    header: &noodles::csi::binning_index::index::Header,
    reference_sequence_name: &str,
) -> io::Result<usize> {
    let Some(id) = header
        .reference_sequence_names()
        .get_index_of(reference_sequence_name.as_bytes())
    else {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "Reference sequence {} not found in index header.",
                reference_sequence_name
            ),
        ));
    };
    Ok(id)
}
