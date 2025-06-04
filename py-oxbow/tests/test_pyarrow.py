import pickle
from typing import Iterable, Iterator

import pyarrow as pa
import pyarrow.dataset as ds

from oxbow._pyarrow import BatchReaderDataset, BatchReaderFragment


def toy_schema():
    return pa.schema(
        [
            ("id", pa.int64()),
            ("name", pa.string()),
            ("age", pa.int32()),
            ("is_active", pa.bool_()),
        ]
    )


def toy_batchreader_builder(columns=None, batch_size=3):
    if columns is None:
        columns = ["id", "name", "age", "is_active"]

    schema = toy_schema()
    schema = pa.schema([schema.field(col) for col in columns])

    data1 = {
        "id": [1, 2, 3],
        "name": ["Alice", "Bob", "Charlie"],
        "age": [25, 30, 35],
        "is_active": [True, False, True],
    }
    arrays1 = {col: pa.array(data1[col]) for col in columns}
    record_batch1 = pa.RecordBatch.from_pydict(arrays1, schema=schema)

    data2 = {
        "id": [4, 5, 6],
        "name": ["David", "Eve", "Frank"],
        "age": [40, 45, 50],
        "is_active": [False, False, True],
    }
    arrays2 = {col: pa.array(data2[col]) for col in columns}
    record_batch2 = pa.RecordBatch.from_pydict(arrays2, schema=schema)

    batch_reader = pa.RecordBatchReader.from_batches(
        schema=schema,
        batches=[record_batch1, record_batch2],
    )
    return batch_reader


def test_fragment():
    fragment = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    assert fragment.schema == toy_schema()
    assert fragment.partition_expression.equals(ds.scalar(True))

    scanner = fragment.scanner()
    assert isinstance(scanner, ds.Scanner)

    batches = fragment.to_batches()
    assert isinstance(batches, Iterable)
    assert isinstance(next(batches), pa.RecordBatch)

    table = fragment.to_table()
    assert isinstance(table, pa.Table)

    head = fragment.head(4)
    assert isinstance(head, pa.Table)

    count = fragment.count_rows()
    assert count == 6


def test_fragment_scanner():
    fragment = BatchReaderFragment(toy_batchreader_builder, toy_schema())

    t1 = fragment.scanner().to_table()
    t2 = fragment.scanner(schema=toy_schema()).to_table()
    assert t1.equals(t2)

    t3 = fragment.scanner(schema=toy_schema(), columns=["id", "name"]).to_table()
    assert len(t3.schema) == 2
    t3 = fragment.scanner(schema=toy_schema(), columns=["id", "name"]).to_table()
    assert len(t3.schema) == 2
    t3 = fragment.scanner(
        schema=toy_schema(), columns=["id", "name"], filter=ds.field("id") > 2
    ).to_table()
    assert len(t3.schema) == 2
    assert len(t3) == 4


def test_fragment_to_batches():
    fragment = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    batches1 = fragment.to_batches()
    batches2 = fragment.to_batches()
    assert isinstance(batches1, Iterator)
    assert isinstance(batches2, Iterator)
    assert isinstance(next(batches1), pa.RecordBatch)
    assert isinstance(next(batches2), pa.RecordBatch)

    batches1 = fragment.to_batches()
    batches2 = fragment.to_batches()
    batches1 = list(batches1)
    batches2 = list(batches2)
    assert batches1 == batches2
    assert len(batches1) == len(batches2) == 2
    assert len(batches1[0].columns) == 4

    batches1 = list(fragment.to_batches(columns=["id", "name"]))
    batches2 = list(fragment.to_batches(columns=["id", "name"]))
    assert batches1 == batches2
    assert len(batches1) == len(batches2) == 2
    assert len(batches1[0].columns) == len(batches1[1].columns) == 2

    batches1 = list(
        fragment.to_batches(columns=["id", "name"], filter=ds.field("id") > 2)
    )
    batches2 = list(
        fragment.to_batches(columns=["id", "name"], filter=ds.field("id") > 2)
    )
    assert batches1 == batches2
    assert len(batches1[0]) == 1
    assert len(batches1[1]) == 3
    assert len(batches1[0].columns) == len(batches1[1].columns) == 2

    batches1 = list(fragment.to_batches(batch_size=3))
    batches2 = list(fragment.to_batches(batch_size=3))
    assert batches1 == batches2
    assert len(batches1) == len(batches2) == 2
    assert len(batches1[0]) == len(batches1[1]) == 3

    batches1 = list(fragment.to_batches(batch_size=42))  # no effect
    batches2 = list(fragment.to_batches(batch_size=42))
    assert batches1 == batches2
    assert len(batches1) == len(batches2) == 2
    assert len(batches1[0]) == len(batches1[1]) == 3


def test_fragment_to_table():
    fragment = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    table1 = fragment.to_table()
    table2 = fragment.to_table()
    assert isinstance(table1, pa.Table)
    assert isinstance(table2, pa.Table)
    assert table1.equals(table2)

    table1 = fragment.to_table(columns=["id", "name"])
    table2 = fragment.to_table(columns=["id", "name"])
    assert table1.equals(table2)
    assert len(table1.schema) == 2
    assert len(table2.schema) == 2

    table1 = fragment.to_table(columns=["id", "name"], filter=ds.field("id") > 2)
    table2 = fragment.to_table(columns=["id", "name"], filter=ds.field("id") > 2)
    assert table1.equals(table2)
    assert len(table1) == 4
    assert len(table2) == 4
    assert len(table1.schema) == 2
    assert len(table2.schema) == 2


def test_fragment_head():
    fragment = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    table1 = fragment.head(4)
    table2 = fragment.head(4)
    assert isinstance(table1, pa.Table)
    assert isinstance(table2, pa.Table)
    assert table1.equals(table2)
    assert len(table1) == 4
    assert len(table2) == 4

    table1 = fragment.head(2)
    table2 = fragment.head(2)
    assert table1.equals(table2)
    assert len(table1) == 2
    assert len(table2) == 2

    table1 = fragment.head(0)
    table2 = fragment.head(0)
    assert table1.equals(table2)
    assert len(table1) == 0
    assert len(table2) == 0


def test_fragment_count_rows():
    fragment = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    count1 = fragment.count_rows()
    count2 = fragment.count_rows()
    assert count1 == count2 == 6

    count1 = fragment.count_rows(filter=ds.field("id") > 2)
    count2 = fragment.count_rows(filter=ds.field("id") > 2)
    assert count1 == count2 == 4

    count1 = fragment.count_rows(filter=ds.field("id") > 6)
    count2 = fragment.count_rows(filter=ds.field("id") > 6)
    assert count1 == count2 == 0


def test_fragment_pickle():
    fragment = BatchReaderFragment(
        toy_batchreader_builder,
        toy_schema(),
        batch_size=3,
    )
    fragment2 = pickle.loads(pickle.dumps(fragment))
    assert isinstance(fragment2, BatchReaderFragment)
    assert fragment2.schema == fragment.schema
    assert fragment2.partition_expression.equals(fragment.partition_expression)
    assert fragment2._batch_size == fragment._batch_size


def test_dataset():
    frag1 = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    frag2 = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    dataset = BatchReaderDataset([frag1, frag2])

    assert dataset.schema == toy_schema()
    assert dataset.partition_expression.equals(ds.scalar(True))

    scanner = dataset.scanner()
    assert isinstance(scanner, ds.Scanner)

    batches = dataset.to_batches()
    assert isinstance(batches, Iterator)

    table = dataset.to_table()
    assert isinstance(table, pa.Table)

    head = dataset.head(4)
    assert isinstance(head, pa.Table)
    assert len(head) == 4

    take = dataset.take([1, 4, 5])
    assert isinstance(take, pa.Table)
    assert len(take) == 3

    count = dataset.count_rows()
    assert count == 12


def test_dataset_filter():
    frag1 = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    frag2 = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    dataset = BatchReaderDataset([frag1, frag2])
    dset = dataset.filter(ds.field("id") > 2)
    assert isinstance(dset, BatchReaderDataset)
    assert dset.count_rows() == 8


def test_dataset_scanner():
    frag1 = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    frag2 = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    dataset = BatchReaderDataset([frag1, frag2])

    t1 = dataset.scanner().to_table()
    t2 = pa.concat_tables([frag1.to_table(), frag2.to_table()])
    assert t1.equals(t2)

    t4 = dataset.scanner(columns=["id", "name"]).to_table()
    assert len(t4.schema) == 2
    assert len(t4) == 12

    t5 = dataset.scanner(columns=["id", "name"], filter=ds.field("id") > 2).to_table()
    assert len(t5.schema) == 2
    assert len(t5) == 8


def test_dataset_to_table():
    frag1 = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    frag2 = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    dataset = BatchReaderDataset([frag1, frag2])

    table1 = dataset.to_table()
    table2 = pa.concat_tables([frag1.to_table(), frag2.to_table()])
    assert table1.equals(table2)
    assert len(table1.schema) == len(table2.schema) == 4

    table1 = dataset.to_table(columns=["id", "name"])
    table2 = pa.concat_tables(
        [frag1.to_table(columns=["id", "name"]), frag2.to_table(columns=["id", "name"])]
    )
    assert table1.equals(table2)
    assert len(table1.schema) == len(table2.schema) == 2

    table1 = dataset.to_table(columns=["id", "name"], filter=ds.field("id") > 2)
    table2 = pa.concat_tables(
        [
            frag1.to_table(columns=["id", "name"], filter=ds.field("id") > 2),
            frag2.to_table(columns=["id", "name"], filter=ds.field("id") > 2),
        ]
    )
    assert table1.equals(table2)
    assert len(table1) == len(table2) == 8
    assert len(table1.schema) == len(table2.schema) == 2


def test_dataset_head():
    frag1 = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    frag2 = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    dataset = BatchReaderDataset([frag1, frag2])

    table1 = dataset.head(0)
    table2 = frag1.head(0)
    assert table1.equals(table2)
    assert len(table1) == 0
    assert len(table2) == 0

    table1 = dataset.head(6)
    table2 = frag1.head(6)
    assert table1.equals(table2)
    assert len(table1) == len(table2) == 6

    table1 = dataset.head(8)
    table2 = pa.concat_tables([frag1.head(6), frag2.head(2)])
    assert table1.equals(table2)
    assert len(table1) == len(table2) == 8


def test_dataset_take():
    frag1 = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    frag2 = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    dataset = BatchReaderDataset([frag1, frag2])

    table = dataset.take([1, 4, 5])
    assert len(table) == 3
    assert table.column("id").to_pylist() == [2, 5, 6]

    table = dataset.take([1, 4, 5, 6, 7, 8, 9])
    assert len(table) == 7
    assert table.column("id").to_pylist() == [2, 5, 6, 1, 2, 3, 4]

    table = dataset.take([1, 4, 5, 6, 7, 8, 9, 10, 11])
    assert len(table) == 9
    assert table.column("id").to_pylist() == [2, 5, 6, 1, 2, 3, 4, 5, 6]


def test_dataset_count_rows():
    frag1 = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    frag2 = BatchReaderFragment(toy_batchreader_builder, toy_schema())
    dataset = BatchReaderDataset([frag1, frag2])

    count1 = dataset.count_rows()
    count2 = frag1.count_rows() + frag2.count_rows()
    assert count1 == count2 == 12

    count1 = dataset.count_rows(filter=ds.field("id") > 2)
    count2 = frag1.count_rows(filter=ds.field("id") > 2) + frag2.count_rows(
        filter=ds.field("id") > 2
    )
    assert count1 == count2 == 8

    count1 = dataset.count_rows(filter=ds.field("id") > 6)
    count2 = frag1.count_rows(filter=ds.field("id") > 6) + frag2.count_rows(
        filter=ds.field("id") > 6
    )
    assert count1 == count2 == 0
