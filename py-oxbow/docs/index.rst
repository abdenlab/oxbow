oxbow
=====

Read specialized bioinformatic file formats as data frames in R, Python, and more.

File formats create a lot of friction for computational biologists. Oxbow is a data unification layer that aims to improve data accessibility and ease of high-performance analytics.

Data I/O is handled in Rust with features exposed to Python and R via Apache Arrow.

Learn more in our recent `blog post <https://open.substack.com/pub/lifeinbytes/p/breaking-out-of-bioinformatic-data-silos?r=jue12&utm_campaign=post&utm_medium=web>`_.

.. grid:: 1 1 2 2
   :gutter: 2

   .. grid-item-card:: Installation
      :link: installation
      :link-type: doc

      Installation instructions for Rust and Python.

   .. grid-item-card:: Quickstart
      :link: quickstart
      :link-type: doc

      Examples of how to use the library's features.

   .. grid-item-card:: API Reference (Python)
      :link: generated/oxbow
      :link-type: doc

      API reference for the Python interface.

   .. grid-item-card:: API Reference (Rust)
      :link: https://docs.rs/oxbow/latest/oxbow/
      :link-type: url

      API reference for the Rust interface.


.. toctree::
   :maxdepth: 2
   :hidden:

   installation

.. toctree::
   :maxdepth: 2
   :hidden:

   quickstart

.. toctree::
   :maxdepth: 3
   :hidden:
   :titlesonly:

   API reference (Python)<generated/oxbow>

.. toctree::
   :maxdepth: 1
   :hidden:

   API reference (Rust)<https://docs.rs/oxbow/latest/oxbow/>