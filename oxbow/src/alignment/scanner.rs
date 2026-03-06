use std::io;

pub mod bam;
pub mod batch_iterator;
pub mod cram;
pub mod sam;

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
