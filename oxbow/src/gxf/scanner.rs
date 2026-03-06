use std::io;

pub mod batch_iterator;
pub mod gff;
pub mod gtf;

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
