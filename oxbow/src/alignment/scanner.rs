use crate::OxbowError;

pub mod bam;
pub mod batch_iterator;
pub mod cram;
pub mod sam;

fn resolve_chrom_id(
    header: &noodles::sam::Header,
    reference_sequence_name: &[u8],
) -> crate::Result<usize> {
    let Some(id) = header
        .reference_sequences()
        .get_index_of(reference_sequence_name)
    else {
        return Err(OxbowError::not_found(format!(
            "Reference sequence {:?} not found in index header.",
            reference_sequence_name
        )));
    };
    Ok(id)
}
