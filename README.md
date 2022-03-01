# Bioframe: Operations on Genomic Interval Dataframes

![Python package](https://github.com/open2c/bioframe/workflows/Python%20package/badge.svg)
[![DOI](https://zenodo.org/badge/69901992.svg)](https://zenodo.org/badge/latestdoi/69901992)
[![Docs status](https://readthedocs.org/projects/bioframe/badge/)](https://bioframe.readthedocs.io/en/latest/)
<img src="./docs/figs/bioframe-logo.png" width=75%> 

Bioframe is a library to enable flexible and scalable operations on genomic interval dataframes in python. Building bioframe directly on top of [pandas](https://pandas.pydata.org/) enables immediate access to a rich set of dataframe operations. Working in python enables rapid visualization (e.g. matplotlib, seaborn) and iteration of genomic analyses.

The philosophy underlying bioframe is to enable flexible operations: instead of creating a function for every possible use-case, we instead encourage users to compose functions to achieve their goals.

Bioframe implements a variety of genomic interval operations directly on dataframes. Bioframe also includes functions for loading diverse genomic data formats, and performing operations on special classes of genomic intervals, including chromosome arms and fixed size bins.

Read the [docs](https://bioframe.readthedocs.io/en/latest/), including the [guide](https://bioframe.readthedocs.io/en/latest/guide-intervalops.html), as well as the [bioframe preprint](https://doi.org/10.1101/2022.02.16.480748) for more information.

If you use ***bioframe*** in your work, please cite:  
*Bioframe: Operations on Genomic Intervals in Pandas Dataframes*. Open2C, Nezar Abdennur, Geoffrey Fudenberg, Ilya Flyamer, Aleksandra A. Galitsyna, Anton Goloborodko, Maxim Imakaev, Sergey V. Venev.
bioRxiv 2022.02.16.480748; doi: https://doi.org/10.1101/2022.02.16.480748


## Installation
The following are required before installing bioframe:
* Python 3.7+
* `numpy`
* `pandas>=1.3`

```sh
pip install bioframe
```

## Interval operations

Key genomic interval operations in bioframe include:
- `closest`: For every interval in a dataframe, find the closest intervals in a second dataframe. 
- `cluster`: Group overlapping intervals in a dataframe into clusters.
- `complement`: Find genomic intervals that are not covered by any interval from a dataframe.
- `overlap`: Find pairs of overlapping genomic intervals between two dataframes. 

Bioframe additionally has functions that are frequently used for genomic interval operations and can be expressed as combinations of these core operations and dataframe operations, including: `coverage`, `expand`, `merge`, `select`, and `subtract`.

To `overlap` two dataframes, call:
```python
import bioframe as bf

bf.overlap(df1, df2)
```

For these two input dataframes, with intervals all on the same chromosome:

<img src="./docs/figs/df1.png" width=60%> 
<img src="./docs/figs/df2.png" width=60%> 

`overlap` will return the following interval pairs as overlaps:

<img src="./docs/figs/overlap_inner_0.png" width=60%> 
<img src="./docs/figs/overlap_inner_1.png" width=60%> 


To `merge` all overlapping intervals in a dataframe, call:
```python
import bioframe as bf

bf.merge(df1)
```

For this input dataframe, with intervals all on the same chromosome:

<img src="./docs/figs/df1.png" width=60%> 

`merge` will return a new dataframe with these merged intervals:

<img src="./docs/figs/merge_df1.png" width=60%> 

See the [guide](https://bioframe.readthedocs.io/en/latest/guide-intervalops.html) for visualizations of other interval operations in bioframe.

## File I/O

Bioframe includes utilities for reading genomic file formats into dataframes and vice versa. One handy function is `read_table` which mirrors pandas’s read_csv/read_table but provides a [`schema`](https://github.com/open2c/bioframe/blob/main/bioframe/io/schemas.py) argument to populate column names for common tabular file formats.

```python
jaspar_url = 'http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/hg38/tsv/MA0139.1.tsv.gz'
ctcf_motif_calls = bioframe.read_table(jaspar_url, schema='jaspar', skiprows=1)
```

## Tutorials
See this [jupyter notebook](https://github.com/open2c/bioframe/tree/master/docs/tutorials/tutorial_assign_motifs_to_peaks.ipynb) for an example of how to assign TF motifs to ChIP-seq peaks using bioframe. 

## Projects currently using bioframe:
* [cooler](https://github.com/open2c/cooler)
* [cooltools](https://github.com/open2c/cooltools)
* yours? :)
