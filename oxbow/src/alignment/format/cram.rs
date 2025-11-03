use std::io::{self, Read, Seek, SeekFrom};
use std::sync::Arc;

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema;
use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use noodles::core::region::Interval;

use crate::alignment::batch_iterator::BatchIterator;
use crate::alignment::model::batch_builder::Push;
use crate::alignment::model::field::DEFAULT_FIELD_NAMES;
use crate::alignment::model::tag::TagScanner;
use crate::alignment::model::BatchBuilder;

/// A CRAM scanner.
///
/// # Examples
///
/// ```no_run
/// use oxbow::alignment::format::cram::Scanner;
/// use std::fs::File;
/// use noodles::fasta::io::indexed_reader::Builder as FastaIndexedReaderBuilder;
/// use noodles::fasta::repository::adapters::IndexedReader as FastaIndexedReaderAdapter;
///
/// let fa_reader = FastaIndexedReaderBuilder::default().build_from_path("reference.fa").unwrap();
/// let adapter = FastaIndexedReaderAdapter::new(fa_reader);
/// let repository = noodles::fasta::Repository::new(adapter);
///
/// let inner = File::open("sample.cram").unwrap();
/// let mut fmt_reader = noodles::cram::io::reader::Builder::default()
///     .set_reference_sequence_repository(repository)
///     .build_from_reader(inner);
/// let header = fmt_reader.read_header().unwrap();
///
/// let scanner = Scanner::new(header);
/// let tag_defs = scanner.tag_defs(&mut fmt_reader, Some(1000)).unwrap();
/// let batches = scanner.scan(fmt_reader, None, Some(tag_defs), None, Some(1000));
/// ```
pub struct Scanner {
    header: noodles::sam::Header,
}

impl Scanner {
    /// Creates a SAM scanner from a SAM header.
    pub fn new(header: noodles::sam::Header) -> Self {
        Self { header }
    }

    /// Returns the SAM header.
    pub fn header(&self) -> noodles::sam::Header {
        self.header.clone()
    }

    /// Returns the reference sequence names.
    pub fn chrom_names(&self) -> Vec<String> {
        self.header
            .reference_sequences()
            .iter()
            .map(|(name, _)| name.to_string())
            .collect()
    }

    /// Returns the reference sequence names and lengths.
    pub fn chrom_sizes(&self) -> Vec<(String, u32)> {
        self.header
            .reference_sequences()
            .iter()
            .map(|(name, r)| (name.to_string(), r.length().get() as u32))
            .collect()
    }

    /// Returns the fixed field names.
    pub fn field_names(&self) -> Vec<String> {
        DEFAULT_FIELD_NAMES.iter().map(|&s| s.to_string()).collect()
    }

    /// Returns the Arrow schema.
    pub fn schema(
        &self,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
    ) -> io::Result<Schema> {
        let header = self.header();
        let batch_builder = BatchBuilder::new(header, fields, tag_defs, 0)?;
        Ok(batch_builder.get_arrow_schema())
    }
}

impl Scanner {
    /// Discovers tag definitions by scanning over records.
    ///
    /// The scan will begin at the current position of the reader and will
    /// move the cursor to the end of the last record scanned.
    pub fn tag_defs<R: Read>(
        &self,
        fmt_reader: &mut noodles::cram::io::Reader<R>,
        scan_rows: Option<usize>,
    ) -> io::Result<Vec<(String, String)>> {
        let header = self.header();
        let records = fmt_reader.records(&header);
        let mut tag_scanner = TagScanner::new();
        match scan_rows {
            None => {
                for result in records {
                    if let Ok(record) = result {
                        tag_scanner.push(&record);
                    } else {
                        eprintln!("Failed to read record");
                    }
                }
            }
            Some(n) => {
                for result in records.take(n) {
                    if let Ok(record) = result {
                        tag_scanner.push(&record);
                    } else {
                        eprintln!("Failed to read record");
                    }
                }
            }
        }
        Ok(tag_scanner.collect())
    }

    /// Returns an iterator over record batches.
    ///
    /// The scan will begin at the current position of the reader and will
    /// move the cursor to the end of the last record scanned.
    pub fn scan<R: Read>(
        &self,
        fmt_reader: noodles::cram::io::Reader<R>,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = BatchBuilder::new(self.header(), fields, tag_defs, batch_size)?;
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }

    /// Returns an iterator over batches of records satisfying a genomic range predicate.
    ///
    /// This operation requires a CRAI Index.
    ///
    /// The scan will traverse one or more CRAM data containers and slices and
    /// filter for records that overlap the given region. The cursor will stop
    /// at the end of the last record scanned.
    #[allow(clippy::too_many_arguments)]
    pub fn scan_query<R: Read + Seek>(
        &self,
        fmt_reader: noodles::cram::io::Reader<R>,
        repo: noodles::fasta::Repository,
        region: noodles::core::Region,
        index: noodles::cram::crai::Index,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let interval = region.interval();

        let batch_builder = BatchBuilder::new(self.header(), fields, tag_defs, batch_size)?;

        let reference_sequence_id = resolve_chrom_id(&self.header, region.name())?;
        let batch_iter = QueryBatchIterator::new(
            fmt_reader,
            self.header(),
            repo,
            index,
            reference_sequence_id,
            interval,
            batch_builder,
            batch_size,
            limit,
        );
        Ok(batch_iter)
    }
}

pub struct QueryBatchIterator<R: std::io::Read + std::io::Seek> {
    query: CramQuery<R>,
    builder: BatchBuilder,
    batch_size: usize,
    limit: usize,
    count: usize,
}

impl<R: std::io::Read + std::io::Seek> QueryBatchIterator<R> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        fmt_reader: noodles::cram::io::Reader<R>,
        header: noodles::sam::Header,
        repo: noodles::fasta::Repository,
        index: noodles::cram::crai::Index,
        reference_sequence_id: usize,
        interval: Interval,
        builder: BatchBuilder,
        batch_size: usize,
        limit: Option<usize>,
    ) -> Self {
        let query = CramQuery::new(
            fmt_reader,
            header,
            index,
            repo,
            reference_sequence_id,
            interval,
        );
        Self {
            query,
            builder,
            batch_size,
            limit: limit.unwrap_or(usize::MAX),
            count: 0,
        }
    }
}

impl<R: std::io::Read + std::io::Seek> RecordBatchReader for QueryBatchIterator<R>
where
    Self: Iterator<Item = Result<RecordBatch, ArrowError>>,
{
    fn schema(&self) -> arrow::datatypes::SchemaRef {
        Arc::new(self.builder.get_arrow_schema())
    }
}

impl<R> Iterator for QueryBatchIterator<R>
where
    R: std::io::Read + std::io::Seek,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut count = 0;
        while count < self.batch_size && self.count < self.limit {
            match self.query.next() {
                Some(Ok(record)) => {
                    match self.builder.push(&record) {
                        Ok(_) => {
                            self.count += 1;
                            count += 1;
                        }
                        Err(e) => return Some(Err(e.into())),
                    };
                }
                Some(Err(e)) => return Some(Err(e.into())),
                None => break,
            }
        }

        if count == 0 {
            None
        } else {
            let batch = self.builder.finish(); // Resets the builder
            Some(batch)
        }
    }
}

pub struct CramQuery<R>
where
    R: Read + Seek,
{
    reader: noodles::cram::io::Reader<R>,
    header: noodles::sam::Header,
    index: std::vec::IntoIter<noodles::cram::crai::Record>,
    repo: noodles::fasta::Repository,
    reference_sequence_id: usize,
    interval: noodles::core::region::Interval,
    records: std::vec::IntoIter<noodles::cram::Record>,
}

impl<R> CramQuery<R>
where
    R: Read + Seek,
{
    pub fn new(
        reader: noodles::cram::io::Reader<R>,
        header: noodles::sam::Header,
        index: Vec<noodles::cram::crai::Record>,
        repo: noodles::fasta::Repository,
        reference_sequence_id: usize,
        interval: noodles::core::region::Interval,
    ) -> Self {
        Self {
            reader,
            header,
            index: index.into_iter(),
            repo,
            reference_sequence_id,
            interval,
            records: Vec::new().into_iter(),
        }
    }

    fn read_next_container(&mut self) -> Option<io::Result<()>> {
        let index_record = self.index.next()?;

        if index_record.reference_sequence_id() != Some(self.reference_sequence_id) {
            return Some(Ok(()));
        }

        if let Err(e) = self.reader.seek(SeekFrom::Start(index_record.offset())) {
            return Some(Err(e));
        }

        let container = match self.reader.read_data_container() {
            Ok(Some(c)) => c,
            Ok(None) => return None,
            Err(e) => return Some(Err(e)),
        };

        let records = container
            .slices()
            .iter()
            .map(|slice| {
                let compression_header = container.compression_header();
                slice.records(compression_header).and_then(|mut records| {
                    slice.resolve_records(
                        &self.repo,
                        &self.header,
                        compression_header,
                        &mut records,
                    )?;

                    Ok(records)
                })
            })
            .collect::<Result<Vec<_>, _>>();

        let records = match records {
            Ok(records) => records,
            Err(e) => return Some(Err(e)),
        };

        self.records = records
            .into_iter()
            .flatten()
            .collect::<Vec<_>>()
            .into_iter();

        Some(Ok(()))
    }
}

impl<R> Iterator for CramQuery<R>
where
    R: Read + Seek,
{
    type Item = io::Result<noodles::cram::Record>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.records.next() {
                Some(r) => {
                    if let (Some(start), Some(end)) = (r.alignment_start(), r.alignment_end()) {
                        let alignment_interval = (start..=end).into();

                        if self.interval.intersects(alignment_interval) {
                            return Some(Ok(r));
                        }
                    }
                }
                None => match self.read_next_container() {
                    Some(Ok(())) => {}
                    Some(Err(e)) => return Some(Err(e)),
                    None => return None,
                },
            }
        }
    }
}

fn resolve_chrom_id(
    header: &noodles::sam::Header,
    reference_sequence_name: &[u8],
) -> io::Result<usize> {
    let Some(id) = header
        .reference_sequences()
        .get_index_of(reference_sequence_name)
    else {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "Reference sequence {:?} not found in index header.",
                reference_sequence_name
            ),
        ));
    };
    Ok(id)
}
