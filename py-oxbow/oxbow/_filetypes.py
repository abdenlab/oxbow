from enum import Enum
from typing import Final

import oxbow.oxbow as ox


class FileType(Enum):
    """
    Enum representing different file types and their associated scanners.
    """

    BAM = ox.PyBamScanner
    BCF = ox.PyBcfScanner
    BED = ox.PyBedScanner
    BigBed = ox.PyBigBedScanner
    BigWig = ox.PyBigWigScanner
    BBIZoom = ox.PyBBIZoomScanner
    FASTA = ox.PyFastaScanner
    FASTQ = ox.PyFastqScanner
    GFF = ox.PyGffScanner
    GTF = ox.PyGtfScanner
    SAM = ox.PySamScanner
    VCF = ox.PyVcfScanner


FILETYPE_BY_NAME: Final[dict[str, FileType]] = {t.name.lower(): t for t in FileType}
