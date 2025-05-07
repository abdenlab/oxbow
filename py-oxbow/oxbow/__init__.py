from oxbow._core.alignment import (
    from_bam as from_bam,
    from_sam as from_sam,
)
from oxbow._core.feature import (
    from_bed as from_bed,
    from_bigbed as from_bigbed,
    from_bigwig as from_bigwig,
    from_gff as from_gff,
    from_gtf as from_gtf,
)
from oxbow._core.sequence import (
    from_fasta as from_fasta,
    from_fastq as from_fastq,
)
from oxbow._core.variant import (
    from_bcf as from_bcf,
    from_vcf as from_vcf,
)

__all__ = [
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
