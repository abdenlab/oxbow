"""
Data file adapters and Arrow scanners.
"""

from oxbow._core.alignment import (
    AlignmentFile as AlignmentFile,
    BamFile as BamFile,
    SamFile as SamFile,
)
from oxbow._core.feature import (
    BedFile as BedFile,
    BigBedFile as BigBedFile,
    BigWigFile as BigWigFile,
    GffFile as GffFile,
    GtfFile as GtfFile,
)
from oxbow._core.sequence import (
    FastaFile as FastaFile,
    FastqFile as FastqFile,
)
from oxbow._core.variant import (
    BcfFile as BcfFile,
    VcfFile as VcfFile,
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

__all__ = [
    "AlignmentFile",
    "BamFile",
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
