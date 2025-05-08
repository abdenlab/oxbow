oxbow
=====

Read specialized bioinformatic file formats as data frames in R, Python, and more.

File formats create a lot of friction for computational biologists. Oxbow is a data unification layer that aims to improve data accessibility and ease of high-performance analytics.

Data I/O is handled in Rust with features exposed to Python and R via Apache Arrow.

Learn more in our recent `blog post <https://open.substack.com/pub/lifeinbytes/p/breaking-out-of-bioinformatic-data-silos?r=jue12&utm_campaign=post&utm_medium=web>`_.

.. grid:: 1 1 2 3
   :gutter: 2

   .. grid-item-card:: How to install
      :link: installation
      :link-type: doc

      Installation instructions for Rust and Python.

   .. grid-item-card:: For users
      :link: userguide
      :link-type: doc

      Examples of how to use the library's features.

   .. grid-item-card:: For developers
      :link: userguide
      :link-type: doc

      How to contribute to Oxbow.

   .. grid-item-card:: API Reference (Python)
      :link: api-python
      :link-type: doc

      API reference for the Python interface.

   .. grid-item-card:: API Reference (Rust)
      :link: https://docs.rs/oxbow/latest/oxbow/
      :link-type: url

      API reference for the Rust interface.

   .. grid-item-card:: API Reference (R)
      :link: api-r
      :link-type: url

      API reference for the R interface.


.. toctree::
   :maxdepth: 2
   :hidden:

   installation

.. toctree::
   :maxdepth: 2
   :hidden:

   userguide

.. toctree::
   :maxdepth: 2
   :hidden:

   contributorguide

.. toctree::
   :maxdepth: 3
   :hidden:

   API reference (Python)<api-python>

.. toctree::
   :maxdepth: 1
   :hidden:

   API reference (Rust)<https://docs.rs/oxbow/latest/oxbow/>

.. toctree::
   :maxdepth: 1
   :hidden:

   API reference (R)<api-r>