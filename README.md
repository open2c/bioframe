# Bioframe: Operations on Genomic Interval Dataframes

<img src="https://github.com/open2c/bioframe/raw/main/docs/figs/bioframe-logo.png" width=75%>

![CI](https://github.com/open2c/bioframe/actions/workflows/ci.yml/badge.svg)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/open2c/bioframe/main.svg)](https://results.pre-commit.ci/latest/github/open2c/bioframe/main)
[![Docs status](https://readthedocs.org/projects/bioframe/badge/)](https://bioframe.readthedocs.io/en/latest/)
[![Paper](https://img.shields.io/badge/DOI-10.1093%2Fbioinformatics%2Fbtae088-blue)](https://doi.org/10.1093/bioinformatics/btae088)
[![Zenodo](https://zenodo.org/badge/69901992.svg)](https://zenodo.org/badge/latestdoi/69901992)
[![Slack](https://img.shields.io/badge/chat-slack-%233F0F3F?logo=slack)](https://bit.ly/open2c-slack)
[![NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](https://www.numfocus.org)

Bioframe enables flexible and scalable operations on genomic interval dataframes in Python.

Bioframe is built directly on top of [Pandas](https://pandas.pydata.org/). Bioframe provides:

* A variety of genomic interval operations that work directly on dataframes.
* Operations for special classes of genomic intervals, including chromosome arms and fixed-size bins.
* Conveniences for diverse tabular genomic data formats and loading genome assembly summary information.

Read the [docs](https://bioframe.readthedocs.io/en/latest/), including the [guide](https://bioframe.readthedocs.io/en/latest/guide-intervalops.html), as well as the [bioframe preprint](https://doi.org/10.1101/2022.02.16.480748) for more information.

Bioframe is an Affiliated Project of [NumFOCUS](https://www.numfocus.org).

## Installation

Bioframe is available on [PyPI](https://pypi.org/project/bioframe/) and [bioconda](https://bioconda.github.io/recipes/bioframe/README.html):

```sh
pip install bioframe
```

## Contributing

Interested in contributing to bioframe? That's great! To get started, check out the [contributing guide](https://github.com/open2c/bioframe/blob/main/CONTRIBUTING.md). Discussions about the project roadmap take place on the [Open2C Slack](https://bit.ly/open2c-slack) and regular developer meetings scheduled there. Anyone can join and participate!


## Interval operations

Key genomic interval operations in bioframe include:
- `overlap`: Find pairs of overlapping genomic intervals between two dataframes.
- `closest`: For every interval in a dataframe, find the closest intervals in a second dataframe.
- `cluster`: Group overlapping intervals in a dataframe into clusters.
- `complement`: Find genomic intervals that are not covered by any interval from a dataframe.

Bioframe additionally has functions that are frequently used for genomic interval operations and can be expressed as combinations of these core operations and dataframe operations, including: `coverage`, `expand`, `merge`, `select`, and `subtract`.

To `overlap` two dataframes, call:
```python
import bioframe as bf

bf.overlap(df1, df2)
```

For these two input dataframes, with intervals all on the same chromosome:

<img src="https://github.com/open2c/bioframe/raw/main/docs/figs/df1.png" width=60%>
<img src="https://github.com/open2c/bioframe/raw/main/docs/figs/df2.png" width=60%>

`overlap` will return the following interval pairs as overlaps:

<img src="https://github.com/open2c/bioframe/raw/main/docs/figs/overlap_inner_0.png" width=60%>
<img src="https://github.com/open2c/bioframe/raw/main/docs/figs/overlap_inner_1.png" width=60%>


To `merge` all overlapping intervals in a dataframe, call:
```python
import bioframe as bf

bf.merge(df1)
```

For this input dataframe, with intervals all on the same chromosome:

<img src="https://github.com/open2c/bioframe/raw/main/docs/figs/df1.png" width=60%>

`merge` will return a new dataframe with these merged intervals:

<img src="https://github.com/open2c/bioframe/raw/main/docs/figs/merge_df1.png" width=60%>

See the [guide](https://bioframe.readthedocs.io/en/latest/guide-intervalops.html) for visualizations of other interval operations in bioframe.

## File I/O

Bioframe includes utilities for reading genomic file formats into dataframes and vice versa. One handy function is `read_table` which mirrors pandas’s read_csv/read_table but provides a [`schema`](https://github.com/open2c/bioframe/blob/main/bioframe/io/schemas.py) argument to populate column names for common tabular file formats.

```python
jaspar_url = 'http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2022/hg38/MA0139.1.tsv.gz'
ctcf_motif_calls = bioframe.read_table(jaspar_url, schema='jaspar', skiprows=1)
```

## Tutorials
See this [jupyter notebook](https://github.com/open2c/bioframe/tree/master/docs/tutorials/tutorial_assign_motifs_to_peaks.ipynb) for an example of how to assign TF motifs to ChIP-seq peaks using bioframe.


## Citing

If you use ***bioframe*** in your work, please cite:

```bibtex
@article{bioframe_2024,
author = {Open2C and Abdennur, Nezar and Fudenberg, Geoffrey and Flyamer, Ilya M and Galitsyna, Aleksandra A and Goloborodko, Anton and Imakaev, Maxim and Venev, Sergey},
doi = {10.1093/bioinformatics/btae088},
journal = {Bioinformatics},
title = {{Bioframe: Operations on Genomic Intervals in Pandas Dataframes}},
year = {2024}
}
```
