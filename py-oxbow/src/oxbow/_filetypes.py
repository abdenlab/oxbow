from enum import Enum
from typing import Final

import oxbow as ox


class FileType(Enum):
    """
    Enum representing different file types and their associated scanners.
    """

    FASTA = ox.PyFastaScanner
    FASTQ = ox.PyFastqScanner
    SAM = ox.PySamScanner
    BAM = ox.PyBamScanner
    VCF = ox.PyVcfScanner
    BCF = ox.PyBcfScanner
    BED = ox.PyBedScanner
    GTF = ox.PyGtfScanner
    GFF = ox.PyGffScanner
    BigWig = ox.PyBigWigScanner
    BigBed = ox.PyBigBedScanner


FILETYPE_BY_NAME: Final[dict[str, FileType]] = {t.name.lower(): t for t in FileType}
