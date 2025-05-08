"""
Data file adapters and Arrow scanners.
"""

from oxbow._core.base import DataSource

from oxbow._core.alignment import (
    AlignmentFile,
    BamFile,
    SamFile,
)
from oxbow._core.feature import (
    BedFile,
    BigBedFile,
    BigWigFile,
    GffFile,
    GtfFile,
)
from oxbow._core.sequence import (
    FastaFile,
    FastqFile,
)
from oxbow._core.variant import (
    BcfFile,
    VcfFile,
)
from oxbow._filetypes import FileType as FileType
from oxbow.oxbow import (
    PyBamScanner,
    PyBcfScanner,
    PyFastaScanner,
    PyFastqScanner,
    PySamScanner,
    PyVcfScanner,
    PyBedScanner,
    PyBigBedScanner,
    PyBigWigScanner,
    PyGffScanner,
    PyGtfScanner,
    PyBBIZoomScanner,
)
from oxbow._filetypes import FileType

__all__ = [
    "AlignmentFile",
    "BamFile",
    "BcfFile",
    "BedFile",
    "BigBedFile",
    "BigWigFile",
    "DataSource",
    "FastaFile",
    "FastqFile",
    "FileType",
    "GffFile",
    "GtfFile",
    "PyBBIZoomScanner",
    "PyBamScanner",
    "PyBcfScanner",
    "PyBedScanner",
    "PyBigBedScanner",
    "PyBigWigScanner",
    "PyFastaScanner",
    "PyFastqScanner",
    "PyGffScanner",
    "PyGtfScanner",
    "PySamScanner",
    "PyVcfScanner",
    "SamFile",
    "VcfFile",
]

AlignmentFile.__module__ = __name__
BamFile.__module__ = __name__
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
