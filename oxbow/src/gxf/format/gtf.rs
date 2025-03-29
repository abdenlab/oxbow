use std::io::{self, BufRead, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema;
use noodles::csi::binning_index;
use noodles::csi::BinningIndex;

use crate::gxf::batch_iterator::{BatchIterator, QueryBatchIterator};
use crate::gxf::model::attribute::AttributeScanner;
use crate::gxf::model::attribute::Push as _;
use crate::gxf::model::field::DEFAULT_FIELD_NAMES;
use crate::gxf::model::BatchBuilder;
use crate::util::query::BgzfChunkReader;

/// A GTF scanner.
///
/// # Examples
///
/// ```no_run
/// use oxbow::gxf::format::gtf::Scanner;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let inner = File::open("sample.gtf").map(BufReader::new).unwrap();
/// let mut fmt_reader = noodles::gtf::io::Reader::new(inner);
///
/// let scanner = Scanner::new(None);
/// let attr_defs = scanner.attribute_defs(&mut fmt_reader, Some(1000)).unwrap();
/// let batches = scanner.scan(fmt_reader, None, Some(attr_defs), None, Some(1000));
/// ```
pub struct Scanner {
    header: Option<binning_index::index::Header>,
}

impl Scanner {
    /// Creates a GTF scanner.
    pub fn new(header: Option<binning_index::index::Header>) -> Self {
        Self { header }
    }

    /// Returns the index header if one was provided.
    pub fn header(&self) -> Option<binning_index::index::Header> {
        self.header.clone()
    }

    /// Returns the reference sequence names if an index header was provided.
    pub fn chrom_names(&self) -> io::Result<Vec<String>> {
        if let Some(header) = &self.header {
            Ok(header
                .reference_sequence_names()
                .iter()
                .map(|name| name.to_string())
                .collect())
        } else {
            Err(io::Error::new(
                io::ErrorKind::NotFound,
                "Index header not found.",
            ))
        }
    }

    /// Returns the fixed field names.
    pub fn field_names(&self) -> Vec<String> {
        DEFAULT_FIELD_NAMES.iter().map(|&s| s.to_string()).collect()
    }

    /// Returns the Arrow schema.
    pub fn schema(
        &self,
        fields: Option<Vec<String>>,
        attr_defs: Option<Vec<(String, String)>>,
    ) -> io::Result<Schema> {
        let batch_builder = BatchBuilder::new(fields, attr_defs, 0)?;
        Ok(batch_builder.get_arrow_schema())
    }
}

impl Scanner {
    /// Discovers attribute definitions by scanning over records.
    ///
    /// The scan will begin at the current position of the reader and will
    /// move the cursor to the end of the last record scanned.
    pub fn attribute_defs<R: BufRead>(
        &self,
        fmt_reader: &mut noodles::gtf::io::Reader<R>,
        scan_rows: Option<usize>,
    ) -> io::Result<Vec<(String, String)>> {
        use noodles::gtf::Line;
        let lines = fmt_reader.lines();
        let mut attr_scanner = AttributeScanner::new();
        match scan_rows {
            None => {
                for line in lines {
                    match line {
                        Ok(line) => {
                            match line {
                                Line::Record(record) => attr_scanner.push(record),
                                Line::Comment(_) => continue,
                            };
                        }
                        Err(e) => eprintln!("Failed to read line: {}", e),
                    }
                }
            }
            Some(n) => {
                for line in lines.take(n) {
                    match line {
                        Ok(line) => {
                            match line {
                                Line::Record(record) => attr_scanner.push(record),
                                Line::Comment(_) => continue,
                            };
                        }
                        Err(e) => eprintln!("Failed to read line: {}", e),
                    }
                }
            }
        }
        Ok(attr_scanner.collect())
    }

    /// Returns an iterator over record batches.
    ///
    /// The scan will begin at the current position of the reader and will
    /// move the cursor to the end of the last record scanned.
    pub fn scan<R: BufRead>(
        &self,
        fmt_reader: noodles::gtf::io::Reader<R>,
        fields: Option<Vec<String>>,
        attr_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = BatchBuilder::new(fields, attr_defs, batch_size)?;
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
    #[allow(clippy::too_many_arguments)]
    pub fn scan_query<R: BufRead + Seek>(
        &self,
        fmt_reader: noodles::gtf::io::Reader<noodles::bgzf::Reader<R>>,
        region: noodles::core::Region,
        index: impl BinningIndex,
        fields: Option<Vec<String>>,
        attr_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let reference_sequence_name = region.name().to_string();
        let interval = region.interval();

        let batch_builder = BatchBuilder::new(fields, attr_defs, batch_size)?;

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
        let fmt_reader = noodles::gtf::io::Reader::new(query_reader);
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
