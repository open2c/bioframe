.. bioframe documentation master file, created by
   sphinx-quickstart on Sat Apr 11 11:44:26 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

bioframe
========

`Bioframe <https://github.com/open2c/bioframe>`_ is a library to enable flexible and scalable operations on genomic interval dataframes in python. Building bioframe directly on top of `pandas <https://pandas.pydata.org/>`_ enables immediate access to a rich set of dataframe operations. Working in python enables rapid visualization and iteration of genomic analyses.


.. toctree::
   :maxdepth: 1
   :caption: Guide

   guide-quickstart
   guide-intervalops.md
   guide-io.ipynb
   guide-performance.ipynb
   guide-recipes.md
   guide-definitions
   guide-specifications
   guide-bedtools

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   tutorials/tutorial_assign_motifs_to_peaks.ipynb
   tutorials/tutorial_assign_peaks_to_genes.ipynb

.. toctree::
   :maxdepth: 3
   :caption: API

   api-construction
   api-validation
   api-intervalops
   api-fileops
   api-resources
   api-extras
   api-vis
   api-lowlevel.md


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
