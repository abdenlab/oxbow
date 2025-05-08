import inspect

import debugpy
import pytest
from wiretap import WireTap


@pytest.fixture
def wiretap():
    def wrapper(target) -> WireTap:
        wiretap = WireTap()
        for cls in inspect.getmro(target):
            for name, _ in inspect.getmembers(
                cls,
                predicate=(
                    lambda x: inspect.isfunction(x)
                    or inspect.ismethod(x)
                    or isinstance(x, property)
                ),
            ):
                if name == "__class__":
                    continue
                else:
                    wiretap.spy(cls, name)
        return wiretap

    yield wrapper


@pytest.hookimpl(tryfirst=True)
def pytest_configure(config):
    """Configures pytest to wait for a debugger to attach if the
    --wait-for-debugger option is set."""
    if config.getoption("--wait-for-debugger"):
        print("Waiting for debugger to attach...")
        debugpy.listen(("0.0.0.0", 5678))
        debugpy.wait_for_client()
        print("Debugger attached.")


def pytest_addoption(parser):
    """Adds custom command-line options to pytest."""
    parser.addoption(
        "-D",
        "--wait-for-debugger",
        action="store_true",
        default=False,
        help="Wait for a debugpy client to attach before running tests",
    )
