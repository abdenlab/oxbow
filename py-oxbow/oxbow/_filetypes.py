from enum import Enum
from typing import Final

import oxbow.oxbow as ox


class FileType(Enum):
    """
    Enum representing different file types and their associated scanners.
    """

    SAM = ox.PySamScanner
    BAM = ox.PyBamScanner


FILETYPE_BY_NAME: Final[dict[str, FileType]] = {t.name.lower(): t for t in FileType}
