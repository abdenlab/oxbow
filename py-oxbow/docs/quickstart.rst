Quickstart
==========

This is a quickstart guide to using Oxbow. It provides several examples of how to use the library's features.

.. code-block:: Python
    
    import oxbow

    # Read a file
    datafile = oxbow.from_bcf("file.bcf", index="file.bcf.csi")

    # Select a dataset
    dataset = datafile.select(samples=["HG00096", "HG00101", "HG00103"], batch_size=10)

    # Convert to a pandas DataFrame
    dataframe = dataset.to_pandas()
