Installation
============

Rust
++++

**From crates.io**

To install the library using the Rust package manager `cargo`, run the following command:

.. code-block:: bash

    cargo install oxbow

**From source**

To install the library from source, clone the repository from GitHub:

.. code-block:: bash

    git clone https://github.com/abdenlab/oxbow.git

The Rust library is located in the `oxbow` directory. To build and install it, ensure you have Rust installed on your system. You can install Rust using `rustup`:

.. code-block:: bash

    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

Once Rust is installed, navigate to the `oxbow` directory and build the library:

.. code-block:: bash

    cd oxbow
    cargo build --release

To install the library system-wide, run:

.. code-block:: bash

    cargo install --path .

Python
++++++

**From PyPI**

To install the library using `pip`, run the following command:

.. code-block:: bash

    pip install oxbow

**From source**

To install the library from source, clone the repository from GitHub:

.. code-block:: bash

    git clone https://github.com/abdenlab/oxbow.git

The Python library is located in the `py-oxbow` directory. To install it, ensure you have Python 3.8 or later and `pip` installed on your system. Then, navigate to the `py-oxbow` directory and install the library:

Using `pip`:

.. code-block:: bash

    cd py-oxbow
    pip install .

Using `uv`:

.. code-block:: bash

    cd py-oxbow
    uv sync
