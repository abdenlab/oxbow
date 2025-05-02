from oxbow._core.alignment import (
    AlignmentFile as AlignmentFile,
    BamFile as BamFile,
    SamFile as SamFile,
    from_bam as from_bam,
    from_sam as from_sam,
)
from oxbow._core.feature import (
    from_bed as from_bed,
    from_bigbed as from_bigbed,
    from_bigwig as from_bigwig,
    from_gff as from_gff,
    from_gtf as from_gtf,
    BedFile as BedFile,
    BigBedFile as BigBedFile,
    BigWigFile as BigWigFile,
    GffFile as GffFile,
    GtfFile as GtfFile,
)
from oxbow._core.sequence import (
    FastaFile as FastaFile,
    FastqFile as FastqFile,
    from_fasta as from_fasta,
    from_fastq as from_fastq,
)
from oxbow._core.variant import (
    BcfFile as BcfFile,
    from_bcf as from_bcf,
    VcfFile as VcfFile,
    from_vcf as from_vcf,
)
from oxbow.oxbow import (
    PyBamScanner as PyBamScanner,
    PyBcfScanner as PyBcfScanner,
    PyFastaScanner as PyFastaScanner,
    PyFastqScanner as PyFastqScanner,
    PySamScanner as PySamScanner,
    PyVcfScanner as PyVcfScanner,
)
from oxbow._filetypes import FileType as FileType
from oxbow._pyarrow import (
    BatchReaderDataset as BatchReaderDataset,
    BatchReaderFragment as BatchReaderFragment,
)

__all__ = [
    "AlignmentFile",
    "BamFile",
    "BatchReaderDataset",
    "BatchReaderFragment",
    "BcfFile",
    "BedFile",
    "BigBedFile",
    "BigWigFile",
    "FileType",
    "FastaFile",
    "FastqFile",
    "GffFile",
    "GtfFile",
    "PyBamScanner",
    "PyBcfScanner",
    "PyFastaScanner",
    "PyFastqScanner",
    "PySamScanner",
    "PyVcfScanner",
    "SamFile",
    "VcfFile",
    "from_bam",
    "from_bcf",
    "from_bed",
    "from_bigbed",
    "from_bigwig",
    "from_fasta",
    "from_fastq",
    "from_gff",
    "from_gtf",
    "from_sam",
    "from_vcf",
]

AlignmentFile.__module__ = __name__
BamFile.__module__ = __name__
BatchReaderDataset.__module__ = __name__
BatchReaderFragment.__module__ = __name__
BcfFile.__module__ = __name__
BedFile.__module__ = __name__
BigBedFile.__module__ = __name__
BigWigFile.__module__ = __name__
FastaFile.__module__ = __name__
FastqFile.__module__ = __name__
FileType.__module__ = __name__
GffFile.__module__ = __name__
GtfFile.__module__ = __name__
PyBamScanner.__module__ = __name__
PyBcfScanner.__module__ = __name__
PyFastaScanner.__module__ = __name__
PyFastqScanner.__module__ = __name__
PySamScanner.__module__ = __name__
PyVcfScanner.__module__ = __name__
SamFile.__module__ = __name__
VcfFile.__module__ = __name__
from_bam.__module__ = __name__
from_bcf.__module__ = __name__
from_bed.__module__ = __name__
from_bigbed.__module__ = __name__
from_bigwig.__module__ = __name__
from_fasta.__module__ = __name__
from_fastq.__module__ = __name__
from_gff.__module__ = __name__
from_gtf.__module__ = __name__
from_sam.__module__ = __name__
from_vcf.__module__ = __name__
