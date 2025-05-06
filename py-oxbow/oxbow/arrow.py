"""
Streamable PyArrow dataset and fragment implementations.
"""

from oxbow._pyarrow import (
    BatchReaderDataset as BatchReaderDataset,
    BatchReaderFragment as BatchReaderFragment,
)

__all__ = [
    "BatchReaderDataset",
    "BatchReaderFragment",
]

BatchReaderDataset.__module__ = __name__
BatchReaderFragment.__module__ = __name__
