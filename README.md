![Bioframe logo](docs/figs/bioframe-logo.png)


Bioframe is a library to enable flexible and scalable operations on genomic interval dataframes in python. Building bioframe directly on top of [pandas](https://pandas.pydata.org/) enables immediate access to a rich set of dataframe operations. Working in python enables rapid visualization (e.g. matplotlib, seaborn) and iteration of genomic analyses.

The philosophy underlying bioframe is to enable flexible operations: instead of creating a function for every possible use-case, we instead encourage users to compose functions to achieve their goals. As a rough rule of thumb, if a function requires three steps and is crucial for genomic interval arithmetic we have included it; conversely if it can be performed in a single line by composing two of the core functions, we have not included it. 

## Core functions
- `closest`: For every interval in a dataframe, find the closest intervals in a second dataframe. 
- `cluster`: Group overlapping intervals in a dataframe into clusters.
- `complement`: Find genomic intervals that are not covered by any interval from a dataframe.
- `overlap`: Find pairs of overlapping genomic intervals between two dataframes. 

Bioframe additionally has functions that are frequently used for genomic interval operations and can be expressed as combinations of these core operations and dataframe operations, including: coverage, expand, merge, select, and subtract.

Bioframe also has functions for loading diverse genomic data formats, and performing operations on special classes of genomic intervals, including chromosome arms and fixed size bins.

Read the [docs](https://bioframe.readthedocs.io/en/genomic_interval_arithmetic/) and explore 

## Genomic interval operations

```python
import bioframe as bf
```

To overlap two dataframes, call
```python
bf.overlap(df1,df2)
```

For these two input dataframes, with intervals all on the same chromosome:


`overlap` will return the following interval pairs:




To merge all overlapping intervals in a dataframe, call:
```python
bf.merge(df1)
```

For this input dataframe, with intervals all on the same chromosome:



`merge` will return a new dataframe with these intervals:


See this [jupyter notebook](https://github.com/mirnylab/bioframe/tree/genomic_interval_arithmetic/docs/notebooks/intervals_tutorials.ipynb) for visualizations of other core bioframe functions.

See this [jupyter notebook](https://github.com/mirnylab/bioframe/tree/genomic_interval_arithmetic/docs/notebooks/tutorial_assign_motifs_to_peaks.ipynb) for an example of how to assign TF motifs to ChIP-seq peaks using bioframe. 


## Requirements
The following are required before installing bioframe:
* Python 3.4+
* `numpy`
* `pandas>=1.0.3`

## Installation
```sh
pip install bioframe
```

## Projects currently using bioframe:
* [cooler](https://github.com/mirnylab/cooler)
* [cooltools](https://github.com/mirnylab/cooltools)
* yours? :)

