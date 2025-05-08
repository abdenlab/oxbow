from __future__ import annotations

from itertools import product
from typing import Generator, Sequence


def swap_quotes(string):
    return string.translate(str.maketrans({"'": '"', '"': "'"}))


class Input:
    """
    Represents a combination of arguments and keyword arguments for
    parameterized tests.

    Examples:
        >>> input = Input(1, 2, a=3, b=4)
        >>> print(input)
        1, 2, a=3, b=4
    """

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def __str__(self):
        return ", ".join(
            (
                *(swap_quotes(repr(v)) for v in self.args),
                *(f"{k}={swap_quotes(repr(v))}" for k, v in self.kwargs.items()),
            )
        )

    @staticmethod
    def permute(*args: Sequence, **kwargs: Sequence) -> Generator[Input]:
        for combination in product(*args, *kwargs.values()):
            arg_combination, kwarg_combination = (
                combination[: len(args)],
                combination[len(args) :],
            )
            yield Input(*arg_combination, **dict(zip(kwargs.keys(), kwarg_combination)))
