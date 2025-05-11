Python API Reference
--------------------

High-level API
=================

Dataset interface
^^^^^^^^^^^^^^^^^^^^

The following functions provide a PyArrow Dataset interface for reading genomic files that may be larger than memory.

.. rubric:: Sequence formats

.. autosummary::
    :toctree: python
    :caption: Sequence formats

    oxbow.from_fasta
    oxbow.from_fastq


.. rubric:: Alignment formats

.. autosummary::
    :toctree: python
    :caption: Alignment formats

    oxbow.from_sam
    oxbow.from_bam


.. rubric:: Variant call formats

.. autosummary::
    :toctree: python
    :caption: Variant call formats

    oxbow.from_vcf
    oxbow.from_bcf


.. rubric:: Interval feature formats

.. autosummary::
    :toctree: python
    :caption: Interval feature formats

    oxbow.from_gtf
    oxbow.from_gff
    oxbow.from_bed


.. rubric:: UCSC BBI formats

.. autosummary::
    :toctree: python
    :caption: UCSC BBI formats

    oxbow.from_bigwig
    oxbow.from_bigbed


Arrow IPC readers
^^^^^^^^^^^^^^^^^^^^

The following functions convert genomic file formats to the Arrow IPC (aka Feather) format as raw bytes. Indexed files support genomic range queries.

.. rubric:: Arrow IPC readers
.. autosummary::
    :toctree: python
    :caption: Arrow IPC readers

    oxbow.read_fasta
    oxbow.read_fastq
    oxbow.read_sam
    oxbow.read_bam
    oxbow.read_vcf
    oxbow.read_bcf
    oxbow.read_gtf
    oxbow.read_gff
    oxbow.read_bed
    oxbow.read_bigwig
    oxbow.read_bigbed

Low-level API
=================

.. rubric:: Scanners

The following classes are wrappers of the Rust "scanner" objects that can read a genomic file format as a stream of Apache Arrow RecordBatches.

.. autosummary::
    :toctree: python
    :caption: Scanners

    oxbow.core.PyFastaScanner
    oxbow.core.PyFastqScanner
    oxbow.core.PySamScanner
    oxbow.core.PyBamScanner
    oxbow.core.PyVcfScanner
    oxbow.core.PyBcfScanner
    oxbow.core.PyGtfScanner
    oxbow.core.PyGffScanner
    oxbow.core.PyBedScanner
    oxbow.core.PyBigWigScanner
    oxbow.core.PyBigBedScanner
    oxbow.core.PyBBIZoomScanner


.. rubric:: PyArrow adapters

The following classes provide a PyArrow Dataset interface over a stream of Arrow RecordBatches supplied by Oxbow's low-level scanners. 
PyArrow Datasets allow working with large datasets that do not fit in memory.

.. autosummary::
    :toctree: python
    :caption: PyArrow adapters

    oxbow.arrow.BatchReaderFragment
    oxbow.arrow.BatchReaderDataset

.. rubric:: Data source classes

.. autosummary::
    :toctree: python
    :caption: DataSources

    oxbow.core.FastaFile
    oxbow.core.FastqFile
    oxbow.core.SamFile
    oxbow.core.BamFile
    oxbow.core.VcfFile
    oxbow.core.BcfFile
    oxbow.core.GtfFile
    oxbow.core.GffFile
    oxbow.core.BedFile
    oxbow.core.BigWigFile
    oxbow.core.BigBedFile