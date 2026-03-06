use std::io::{self, BufRead, Read, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::{Schema, SchemaRef};
use noodles::bgzf::VirtualPosition;
use noodles::csi::binning_index;
use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles::csi::BinningIndex;

use crate::batch::RecordBatchBuilder as _;
use crate::gxf::model::attribute::AttributeScanner;
use crate::gxf::model::attribute::Push as _;
use crate::gxf::model::field::DEFAULT_FIELD_NAMES;
use crate::gxf::model::BatchBuilder;
use crate::gxf::scanner::batch_iterator::{BatchIterator, QueryBatchIterator};
use crate::util::query::{BgzfChunkReader, ByteRangeReader};

/// A GFF scanner.
///
/// Schema parameters (fields, attribute definitions) are declared at
/// construction time. Scan methods accept only column projection,
/// batch size, and limit.
///
/// # Examples
///
/// ```no_run
/// use oxbow::gxf::scanner::gff::Scanner;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let inner = File::open("sample.gff").map(BufReader::new).unwrap();
/// let mut fmt_reader = noodles::gff::io::Reader::new(inner);
///
/// let attr_defs = Scanner::attribute_defs(&mut fmt_reader, Some(1000)).unwrap();
/// let scanner = Scanner::new(None, None, Some(attr_defs)).unwrap();
/// let batches = scanner.scan(fmt_reader, None, None, Some(1000));
/// ```
pub struct Scanner {
    header: Option<binning_index::index::Header>,
    fields: Option<Vec<String>>,
    attr_defs: Option<Vec<(String, String)>>,
    schema: SchemaRef,
}

impl Scanner {
    /// Creates a GFF scanner from schema parameters.
    ///
    /// The schema is validated and cached at construction time.
    pub fn new(
        header: Option<binning_index::index::Header>,
        fields: Option<Vec<String>>,
        attr_defs: Option<Vec<(String, String)>>,
    ) -> io::Result<Self> {
        let batch_builder = BatchBuilder::new(fields.clone(), attr_defs.clone(), 0)?;
        let schema = batch_builder.schema();
        Ok(Self {
            header,
            fields,
            attr_defs,
            schema,
        })
    }

    /// Returns the index header if one was provided.
    pub fn header(&self) -> Option<&binning_index::index::Header> {
        self.header.as_ref()
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
    pub fn schema(&self) -> &Schema {
        &self.schema
    }

    /// Builds a BatchBuilder applying column projection.
    ///
    /// Returns an error if any requested column is not in the declared schema.
    fn build_batch_builder(
        &self,
        columns: Option<Vec<String>>,
        capacity: usize,
    ) -> io::Result<BatchBuilder> {
        match columns {
            None => BatchBuilder::new(self.fields.clone(), self.attr_defs.clone(), capacity),
            Some(cols) => {
                let schema_names: Vec<&str> = self
                    .schema
                    .fields()
                    .iter()
                    .map(|f| f.name().as_str())
                    .collect();

                let unknown: Vec<&str> = cols
                    .iter()
                    .filter(|c| !schema_names.iter().any(|s| s.eq_ignore_ascii_case(c)))
                    .map(|c| c.as_str())
                    .collect();
                if !unknown.is_empty() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "Unknown columns: {:?}. Available columns: {:?}",
                            unknown, schema_names
                        ),
                    ));
                }

                let declared_field_names: Vec<String> = self.fields.clone().unwrap_or_else(|| {
                    DEFAULT_FIELD_NAMES.iter().map(|&s| s.to_string()).collect()
                });
                let projected_fields: Vec<String> = declared_field_names
                    .into_iter()
                    .filter(|name| cols.iter().any(|c| c.eq_ignore_ascii_case(name)))
                    .collect();

                let attr_defs = if cols.iter().any(|c| c == "attributes") {
                    self.attr_defs.clone()
                } else {
                    None
                };

                BatchBuilder::new(Some(projected_fields), attr_defs, capacity)
            }
        }
    }
}

impl Scanner {
    /// Discovers attribute definitions by scanning over records.
    ///
    /// The scan will begin at the current position of the reader and will
    /// move the cursor to the end of the last record scanned.
    pub fn attribute_defs<R: BufRead>(
        fmt_reader: &mut noodles::gff::io::Reader<R>,
        scan_rows: Option<usize>,
    ) -> io::Result<Vec<(String, String)>> {
        let lines = fmt_reader.lines();
        let mut attr_scanner = AttributeScanner::new();
        match scan_rows {
            None => {
                for line in lines {
                    match line {
                        Ok(line) => match line.as_record() {
                            Some(result) => attr_scanner.push(result?),
                            None => continue,
                        },
                        Err(e) => eprintln!("Failed to read line: {}", e),
                    }
                }
            }
            Some(n) => {
                for line in lines.take(n) {
                    match line {
                        Ok(line) => match line.as_record() {
                            Some(result) => attr_scanner.push(result?),
                            None => continue,
                        },
                        Err(e) => eprintln!("Failed to read line: {}", e),
                    }
                }
            }
        }
        Ok(attr_scanner.collect())
    }

    /// Returns an iterator yielding record batches.
    pub fn scan<R: BufRead>(
        &self,
        fmt_reader: noodles::gff::io::Reader<R>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }

    /// Returns an iterator yielding record batches satisfying a genomic range query.
    pub fn scan_query<R: BufRead + Seek>(
        &self,
        fmt_reader: noodles::gff::io::Reader<noodles::bgzf::io::Reader<R>>,
        region: noodles::core::Region,
        index: impl BinningIndex,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let reference_sequence_name = region.name().to_string();
        let interval = region.interval();

        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        let Some(header) = index.header() else {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                "Index header not found.",
            ));
        };
        let reference_sequence_id = super::resolve_chrom_id(header, &reference_sequence_name)?;
        let chunks = index.query(reference_sequence_id, interval)?;
        let bgzf_reader = fmt_reader.into_inner();
        let query_reader = BgzfChunkReader::new(bgzf_reader, chunks);
        let fmt_reader = noodles::gff::io::Reader::new(query_reader);
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

    /// Returns an iterator yielding record batches from specified byte ranges.
    pub fn scan_byte_ranges<R: BufRead + Seek>(
        &self,
        fmt_reader: noodles::gff::io::Reader<R>,
        byte_ranges: Vec<(u64, u64)>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        let inner_reader = fmt_reader.into_inner();
        let range_reader = ByteRangeReader::new(inner_reader, byte_ranges);
        let fmt_reader = noodles::gff::io::Reader::new(range_reader);
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }

    /// Returns an iterator yielding record batches from specified virtual position ranges.
    pub fn scan_virtual_ranges<R: Read + Seek>(
        &self,
        fmt_reader: noodles::gff::io::Reader<noodles::bgzf::io::Reader<R>>,
        vpos_ranges: Vec<(VirtualPosition, VirtualPosition)>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        let chunks: Vec<Chunk> = vpos_ranges
            .into_iter()
            .map(|(start, end)| Chunk::new(start, end))
            .collect();

        let bgzf_reader = fmt_reader.into_inner();
        let range_reader = BgzfChunkReader::new(bgzf_reader, chunks);
        let fmt_reader = noodles::gff::io::Reader::new(range_reader);
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }
}
