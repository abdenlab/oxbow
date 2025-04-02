from typing import List, Literal, Tuple

class PyBamScanner:
    def __init__(
        self, obj: object, compressed: bool | None = True
    ) -> None: ...
    def chrom_names(self) -> List[str]: ...
    def chrom_sizes(self) -> List[Tuple[str, int]]: ...
    def field_names(self) -> List[str]: ...
    def tag_defs(
        self, scan_rows: int | None = None
    ) -> List[Tuple[str, str]]: ...
    def schema(
        self,
        fields: List[str] | None = None,
        tag_defs: List[Tuple[str, str]] | None = None,
    ) -> object: ...
    def scan(
        self,
        n_records: int | None = None,
        fields: List[str] | None = None,
        tag_defs: List[Tuple[str, str]] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_until(
        self,
        position: int | None = None,
        fields: List[str] | None = None,
        tag_defs: List[Tuple[str, str]] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_until_vpos(
        self,
        position: Tuple[int, int],
        fields: List[str] | None = None,
        tag_defs: List[Tuple[str, str]] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_query(
        self,
        region: str,
        index: object | None = None,
        fields: List[str] | None = None,
        tag_defs: List[Tuple[str, str]] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_unmapped(
        self,
        index: object | None = None,
        fields: List[str] | None = None,
        tag_defs: List[Tuple[str, str]] | None = None,
        batch_size: int | None = None,
    ) -> object: ...

class PyBBIZoomScanner:
    def __init__(
        self, reader: object, bbi_type: object, zoom_level: int
    ) -> None: ...
    def field_names(self) -> List[str]: ...
    def schema(self, fields: List[str] | None = None) -> object: ...
    def scan(
        self,
        n_records: int | None = None,
        fields: List[str] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_query(
        self,
        region: str,
        fields: List[str] | None = None,
        batch_size: int | None = None,
    ) -> object: ...

class PyBedScanner:
    def __init__(
        self, obj: object, bed_schema: str, compressed: bool | None = False
    ) -> None: ...
    def field_names(self) -> List[str]: ...
    def schema(self, fields: List[str] | None = None) -> object: ...
    def scan(
        self,
        n_records: int | None = None,
        fields: List[str] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_until(
        self,
        position: int | None = None,
        fields: List[str] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_until_vpos(
        self,
        position: Tuple[int, int],
        fields: List[str] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_query(
        self,
        region: str,
        index: object | None = None,
        fields: List[str] | None = None,
        batch_size: int | None = None,
    ) -> object: ...

class PyBcfScanner:
    def __init__(
        self, obj: object, compressed: bool | None = True
    ) -> None: ...
    def chrom_names(self) -> List[str]: ...
    def chrom_sizes(self) -> List[Tuple[str, int]]: ...
    def field_names(self) -> List[str]: ...
    def info_fields(self) -> List[Tuple[str, str, str]]: ...
    def genotype_fields(self) -> List[Tuple[str, str, str]]: ...
    def samples(self) -> List[str]: ...
    def schema(
        self,
        fields: List[str] | None = None,
        info_fields: List[Tuple[str, str, str]] | None = None,
        genotype_fields: List[Tuple[str, str, str]] | None = None,
        samples: List[str] | None = None,
        genotype_by: str | None = None,
    ) -> object: ...
    def scan(
        self,
        n_records: int | None = None,
        fields: List[str] | None = None,
        info_fields: List[Tuple[str, str, str]] | None = None,
        genotype_fields: List[Tuple[str, str, str]] | None = None,
        samples: List[str] | None = None,
        genotype_by: str | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_until(
        self,
        position: int | None = None,
        fields: List[str] | None = None,
        info_fields: List[Tuple[str, str, str]] | None = None,
        genotype_fields: List[Tuple[str, str, str]] | None = None,
        samples: List[str] | None = None,
        genotype_by: str | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_until_vpos(
        self,
        position: Tuple[int, int],
        fields: List[str] | None = None,
        info_fields: List[Tuple[str, str, str]] | None = None,
        genotype_fields: List[Tuple[str, str, str]] | None = None,
        samples: List[str] | None = None,
        genotype_by: str | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_query(
        self,
        region: str,
        index: object | None = None,
        fields: List[str] | None = None,
        info_fields: List[Tuple[str, str, str]] | None = None,
        genotype_fields: List[Tuple[str, str, str]] | None = None,
        samples: List[str] | None = None,
        genotype_by: str | None = None,
        batch_size: int | None = None,
    ) -> object: ...

class PyBigBedScanner:
    def __init__(self, obj: object, schema: Literal["autosql", "bedgraph"] | str | None = None) -> None: ...
    def chrom_names(self) -> List[str]: ...
    def chrom_sizes(self) -> List[Tuple[str, int]]: ...
    def field_names(self) -> List[str]: ...
    def schema(self, fields: List[str] | None = None) -> object: ...
    def read_autosql(self) -> str: ...
    def zoom_levels(self) -> List[int]: ...
    def get_zoom(self, zoom_level: int) -> "PyBBIZoomScanner": ...
    def scan(
        self,
        n_records: int | None = None,
        fields: List[str] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_query(
        self,
        region: str,
        fields: List[str] | None = None,
        batch_size: int | None = None,
    ) -> object: ...

class PyBigWigScanner:
    def __init__(self, obj: object) -> None: ...
    def chrom_names(self) -> List[str]: ...
    def chrom_sizes(self) -> List[Tuple[str, int]]: ...
    def field_names(self) -> List[str]: ...
    def schema(self, fields: List[str] | None = None) -> object: ...
    def zoom_levels(self) -> List[int]: ...
    def get_zoom(self, zoom_level: int) -> "PyBBIZoomScanner": ...
    def scan(
        self,
        n_records: int | None = None,
        fields: List[str] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_query(
        self,
        region: str,
        fields: List[str] | None = None,
        batch_size: int | None = None,
    ) -> object: ...

class PyFastqScanner:
    def __init__(
        self, obj: object, compressed: bool | None = False
    ) -> None: ...
    def field_names(self) -> List[str]: ...
    def schema(self, fields: List[str] | None = None) -> object: ...
    def scan(
        self,
        n_records: int | None = None,
        fields: List[str] | None = None,
        batch_size: int | None = None,
    ) -> object: ...

class PyFastaScanner:
    def __init__(
        self, obj: object, compressed: bool | None = False
    ) -> None: ...
    def field_names(self) -> List[str]: ...
    def schema(self, fields: List[str] | None = None) -> object: ...
    def scan(
        self,
        n_records: int | None = None,
        fields: List[str] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_query(
        self,
        regions: List[str],
        index: object | None = None,
        gzi: object | None = None,
        fields: List[str] | None = None,
        batch_size: int | None = None,
    ) -> object: ...

class PyGffScanner:
    def __init__(
        self, obj: object, compressed: bool | None = False
    ) -> None: ...
    def field_names(self) -> List[str]: ...
    def attribute_defs(
        self, scan_rows: int | None = None
    ) -> List[Tuple[str, str]]: ...
    def schema(
        self,
        fields: List[str] | None = None,
        attribute_defs: List[Tuple[str, str]] | None = None,
    ) -> object: ...
    def scan(
        self,
        n_records: int | None = None,
        fields: List[str] | None = None,
        attribute_defs: List[Tuple[str, str]] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_until(
        self,
        position: int | None = None,
        fields: List[str] | None = None,
        attribute_defs: List[Tuple[str, str]] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_until_vpos(
        self,
        position: Tuple[int, int],
        fields: List[str] | None = None,
        attribute_defs: List[Tuple[str, str]] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_query(
        self,
        region: str,
        index: object | None = None,
        fields: List[str] | None = None,
        attribute_defs: List[Tuple[str, str]] | None = None,
        batch_size: int | None = None,
    ) -> object: ...

class PyGtfScanner:
    def __init__(
        self, obj: object, compressed: bool | None = False
    ) -> None: ...
    def field_names(self) -> List[str]: ...
    def attribute_defs(
        self, scan_rows: int | None = None
    ) -> List[Tuple[str, str]]: ...
    def schema(
        self,
        fields: List[str] | None = None,
        attribute_defs: List[Tuple[str, str]] | None = None,
    ) -> object: ...
    def scan(
        self,
        n_records: int | None = None,
        fields: List[str] | None = None,
        attribute_defs: List[Tuple[str, str]] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_until(
        self,
        position: int | None = None,
        fields: List[str] | None = None,
        attribute_defs: List[Tuple[str, str]] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_until_vpos(
        self,
        position: Tuple[int, int],
        fields: List[str] | None = None,
        attribute_defs: List[Tuple[str, str]] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_query(
        self,
        region: str,
        index: object | None = None,
        fields: List[str] | None = None,
        attribute_defs: List[Tuple[str, str]] | None = None,
        batch_size: int | None = None,
    ) -> object: ...

class PySamScanner:
    def __init__(
        self, obj: object, compressed: bool | None = False
    ) -> None: ...
    def chrom_names(self) -> List[str]: ...
    def chrom_sizes(self) -> List[Tuple[str, int]]: ...
    def field_names(self) -> List[str]: ...
    def tag_defs(
        self, scan_rows: int | None = None
    ) -> List[Tuple[str, str]]: ...
    def schema(
        self,
        fields: List[str] | None = None,
        tag_defs: List[Tuple[str, str]] | None = None,
    ) -> object: ...
    def scan(
        self,
        n_records: int | None = None,
        fields: List[str] | None = None,
        tag_defs: List[Tuple[str, str]] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_until(
        self,
        position: int | None = None,
        fields: List[str] | None = None,
        tag_defs: List[Tuple[str, str]] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_until_vpos(
        self,
        position: Tuple[int, int],
        fields: List[str] | None = None,
        tag_defs: List[Tuple[str, str]] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_query(
        self,
        region: str,
        index: object | None = None,
        fields: List[str] | None = None,
        tag_defs: List[Tuple[str, str]] | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_unmapped(
        self,
        index: object | None = None,
        fields: List[str] | None = None,
        tag_defs: List[Tuple[str, str]] | None = None,
        batch_size: int | None = None,
    ) -> object: ...

class PyVcfScanner:
    def __init__(
        self, obj: object, compressed: bool | None = False
    ) -> None: ...
    def chrom_names(self) -> List[str]: ...
    def chrom_sizes(self) -> List[Tuple[str, int]]: ...
    def field_names(self) -> List[str]: ...
    def info_fields(self) -> List[Tuple[str, str, str]]: ...
    def genotype_fields(self) -> List[Tuple[str, str, str]]: ...
    def samples(self) -> List[str]: ...
    def schema(
        self,
        fields: List[str] | None = None,
        info_fields: List[Tuple[str, str, str]] | None = None,
        genotype_fields: List[Tuple[str, str, str]] | None = None,
        samples: List[str] | None = None,
        genotype_by: str | None = None,
    ) -> object: ...
    def scan(
        self,
        n_records: int | None = None,
        fields: List[str] | None = None,
        info_fields: List[Tuple[str, str, str]] | None = None,
        genotype_fields: List[Tuple[str, str, str]] | None = None,
        samples: List[str] | None = None,
        genotype_by: str | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_until(
        self,
        position: int | None = None,
        fields: List[str] | None = None,
        info_fields: List[Tuple[str, str, str]] | None = None,
        genotype_fields: List[Tuple[str, str, str]] | None = None,
        samples: List[str] | None = None,
        genotype_by: str | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_until_vpos(
        self,
        position: Tuple[int, int],
        fields: List[str] | None = None,
        info_fields: List[Tuple[str, str, str]] | None = None,
        genotype_fields: List[Tuple[str, str, str]] | None = None,
        samples: List[str] | None = None,
        genotype_by: str | None = None,
        batch_size: int | None = None,
    ) -> object: ...
    def scan_query(
        self,
        region: str,
        index: object | None = None,
        fields: List[str] | None = None,
        info_fields: List[Tuple[str, str, str]] | None = None,
        genotype_fields: List[Tuple[str, str, str]] | None = None,
        samples: List[str] | None = None,
        genotype_by: str | None = None,
        batch_size: int | None = None,
    ) -> object: ...
