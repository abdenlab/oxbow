from importlib.metadata import PackageNotFoundError, version

from oxbow._core.alignment import (
    from_bam,
    from_sam,
)
from oxbow._core.bbi import (
    from_bigbed,
    from_bigwig,
)
from oxbow._core.bed import (
    from_bed,
)
from oxbow._core.gxf import (
    from_gff,
    from_gtf,
)
from oxbow._core.sequence import (
    from_fasta,
    from_fastq,
)
from oxbow._core.variant import (
    from_bcf,
    from_vcf,
)
from oxbow.oxbow import (
    read_bam,
    read_bcf,
    read_bed,
    read_bigbed,
    read_bigwig,
    read_fasta,
    read_fastq,
    read_gff,
    read_gtf,
    read_sam,
    read_vcf,
)

try:
    __version__ = version("oxbow")
except PackageNotFoundError:
    __version__ = "unknown"

__all__ = [
    "__version__",
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
    "read_fasta",
    "read_fastq",
    "read_sam",
    "read_bam",
    "read_bcf",
    "read_vcf",
    "read_bed",
    "read_bigbed",
    "read_bigwig",
    "read_gff",
    "read_gtf",
]

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

del version, PackageNotFoundError
