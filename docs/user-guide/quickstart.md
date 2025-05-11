# Quickstart

This is a quickstart guide to using Oxbow. It provides several examples of how to use the library's features.

## Reading a File

Files can be read a number of ways. The simplest way is to use the convenience function associated with your file type. The returned `DataSource` object can be used to access the data in the file.

Below are examples of reading a BCF file using the `from_bcf` function, the `DataSource.from_bcf` class method, and the `VariantFile.from_bcf` class method. You can also instantiate a `BcfFile` object directly.

```python    
import oxbow

# Using the `from_bcf` function to read a BCF file
datafile = oxbow.from_bcf("file.bcf", index="file.bcf.csi")

# Using the `DataSource.from_bcf` class method
datafile = oxbow.DataSource.from_bcf("file.bcf", index="file.bcf.csi")

# Using the `VariantFile.from_bcf` class method
datafile = oxbow.VariantFile.from_bcf("file.bcf", index="file.bcf.csi")

# Instantiating a `BcfFile` object directly
datafile = oxbow.BcfFile("file.bcf", index="file.bcf.csi")
```
