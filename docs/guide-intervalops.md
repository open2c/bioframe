---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.3
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Genomic interval operations 

```{eval-rst}
This guide provides an introdution into how to use bioframe to perform genomic interval operations. For the full list of genomic interval operations, see the :ref:`API reference<API_ops>`.

The following modules are used in this guide:
```
```{code-cell} ipython3
import itertools

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

import bioframe as bf
import bioframe.vis
```

### DataFrames & BedFrames

```{eval-rst}
The core objects in bioframe are pandas DatFrames of genomic intervals, or BedFrames. These can either be defined directly with :py:class:`pandas.DataFrame`:
```
```{code-cell} ipython3
df1 = pd.DataFrame([
    ['chr1', 1, 5],
    ['chr1', 3, 8],
    ['chr1', 8, 10],
    ['chr1', 12, 14]],
    columns=['chrom', 'start', 'end']
)
```
```{eval-rst}
Or via functions in :mod:`bioframe.core.construction`, e.g.:
```
```{code-cell} ipython3
df2 = bioframe.from_list(
    [['chr1', 4, 8],
     ['chr1', 10, 11]], 
    name_col='chrom')
```
```{eval-rst}
Or ingested from datasets and databases with functions in :mod:`bioframe.io.fileops`
and :mod:`bioframe.io.resources`.
```

```{eval-rst}
BedFrames satisfy the following properties:  

- chrom, start, end columns  
- columns have valid dtypes (object/string/categorical, int, int)  
- all starts < ends.  

Whether a dataframe satisfies these properties can be checked with :func:`bioframe.core.checks.is_bedframe`:
```
```{code-cell} ipython3
bioframe.is_bedframe(df2)
```
```{eval-rst}
See :mod:`bioframe.core.checks` for other functions that test properties of BedFrames and :ref:`Technical Notes<Technical_Notes>` for detailed definitions.
```

```{eval-rst}
:py:mod:`bioframe.vis` provides plotting utilities for intervals:
```
```{code-cell} ipython3
bf.vis.plot_intervals(df1, show_coords=True, xlim=(0,16))
plt.title('bedFrame1 intervals');

bf.vis.plot_intervals(df2, show_coords=True, xlim=(0,16), colors='lightpink')
plt.title('bedFrame2 intervals');
```

### Overlap
```{eval-rst}
Calculating the overlap between two sets of genomic intervals is a crucial genomic interval operation.

Using :func:`bioframe.overlap`, we can see the two dataframes defined above, ``df1`` and ``df2``, contain two pairs of overlapping intervals:
```
```{code-cell} ipython3
overlapping_intervals = bf.overlap(df1, df2, how='inner', suffixes=('_1','_2'))
display(overlapping_intervals)
```
```{code-cell} ipython3
for i, reg_pair in overlapping_intervals.iterrows(): 
    bf.vis.plot_intervals_arr(
        starts = [reg_pair.start_1,reg_pair.start_2],
        ends = [reg_pair.end_1,reg_pair.end_2],
        colors = ['skyblue', 'lightpink'],
        levels = [2,1],
        xlim = (0,16),
        show_coords = True)
    plt.title(f'overlapping pair #{i}')
```
Note that we passed custom suffixes for the outputs (defaults are ``suffixes=("","_")``),
as well as a custom overlap mode (``how='inner'``). The default overlap mode, ``how='left'`` returns each interval in ``df1`` whether or not it overlaps an interval in ``df2``.
```{code-cell} ipython3
overlapping_intervals = bf.overlap(df1, df2)
display(overlapping_intervals)
```


### Cluster
```{eval-rst}
It is often useful to find overlapping intervals within a single set of genomic intervals. In `bioframe`, this is achieved with :func:`bioframe.cluster`. This function returns a DataFrame where subsets of overlapping intervals are assigned to the same group, reported in a new column.

To demonstrate the usage of :func:`bioframe.cluster`, we use the same ``df1`` as above:
```
```{code-cell} ipython3
df1 = pd.DataFrame([
    ['chr1', 1, 5],
    ['chr1', 3, 8],
    ['chr1', 8, 10],
    ['chr1', 12, 14],
    ],
    columns=['chrom', 'start', 'end']
)

bf.vis.plot_intervals(df1, show_coords=True, xlim=(0,16))
```

Cluster returns a DataFrame where each interval is assigned to a group:
```{code-cell} ipython3
df_annotated = bf.cluster(df1, min_dist=0)
display(df_annotated)
bf.vis.plot_intervals(df_annotated, labels=df_annotated['cluster'], show_coords=True, xlim=(0,16))
```

Note that using ``min_dist=0`` and ``min_dist=None`` give different results, as the latter only clusters overlapping intervals and not adjacent intervals:
```{code-cell} ipython3
df_annotated = bf.cluster(df1, min_dist=None)
display(df_annotated)
bf.vis.plot_intervals(df_annotated, labels=df_annotated['cluster'], show_coords=True, xlim=(0,16))
```

Extending the minimum distance to two (``min_dist=2``) makes all intervals part of the same cluster "0":
```{code-cell} ipython3
df_annotated = bf.cluster(df1, min_dist=2)
display(df_annotated)
bf.vis.plot_intervals(df_annotated, labels=df_annotated['cluster'], show_coords=True, xlim=(0,16))
```

### Merge
```{eval-rst}
Instead of returning cluster assignments, :func:`bioframe.merge` returns a new dataframe of merged genomic intervals. As with :func:`bioframe.cluster`, using ``min_dist=0`` and ``min_dist=None`` gives different results. 

If ``min_dist=0``, this returns a dataframe of two intervals:
```
```{code-cell} ipython3
df_merged = bf.merge(df1, min_dist=0)

display(df_merged)
bf.vis.plot_intervals(df_merged, show_coords=True, xlim=(0,16))
```

If ``min_dist=None``, this returns a dataframe of three intervals:
```{code-cell} ipython3
df_merged = bf.merge(df1, min_dist=None)
display(df_merged)
bf.vis.plot_intervals(df_merged, show_coords=True, xlim=(0,16))
```

### Closest
```{eval-rst}
In genomics, it is often useful not only to find features that overlap, but also features that are nearby along the genome. In bioframe, this is achieved using :func:`bioframe.closest`.
```
```{code-cell} ipython3
closest_intervals = bf.closest(df1, df2, suffixes=('_1','_2'))
display(closest_intervals)
```
```{code-cell} ipython3
for i, reg_pair in closest_intervals.iterrows(): 
    bf.vis.plot_intervals_arr(
        starts = [reg_pair.start_1,reg_pair.start_2],
        ends = [reg_pair.end_1,reg_pair.end_2],
        colors = ['skyblue', 'lightpink'],
        levels = [2,1],
        xlim = (0,16),
        show_coords = True)
    plt.title(f'closest pair #{i}')
```

```{eval-rst}
By default, :func:`bioframe.closest` reports overlapping intervals. This behavior can be modified, however, by passing ``ignore_overlap=True``. Note the closest pair #2 and #3, which did not overlap, remain the same: 
```
```{code-cell} ipython3
closest_intervals = bf.closest(df1, df2, ignore_overlaps=True, suffixes=('_1','_2'))
for i, reg_pair in closest_intervals.iterrows(): 
    bf.vis.plot_intervals_arr(
        starts = [reg_pair.start_1,reg_pair.start_2],
        ends = [reg_pair.end_1,reg_pair.end_2],
        colors = ['skyblue', 'lightpink'],
        levels = [2,1],
        xlim = (0,16),
        show_coords = True)
    plt.title(f'closest pair #{i}')
```

```{eval-rst}
Closest intervals within a single DataFrame can be found simply by passing a single dataframe to :func:`bioframe.closest`. The number of closest intervals to report per query interval can be adjusted with ``k``. 
```
```{code-cell} ipython3
bf.closest(df1, k=2)
```

### Coverage & Count Overlaps
```{eval-rst}
For two sets of genomic features, it is often useful to calculate the number of basepairs covered and the number of overlapping intervals. While these are fairly straightforward to compute from the output of :func:`bioframe.overlap` with pandas.groupby and column renaming, since these are very frequently used, they are provided as core bioframe functions.
```

```{code-cell} ipython3
df1_coverage = bf.coverage(df1, df2)
display(df1_coverage)
```
```{code-cell} ipython3
df1_coverage_and_count = bf.count_overlaps(df1_coverage, df2)
display(df1_coverage_and_count)
```

This plot shows the coverage and number of overlaps for intervals in ``df1`` by ``df2``:
```{code-cell} ipython3
bf.vis.plot_intervals(
    df1_coverage_and_count, 
    show_coords=True, xlim=(0,16), 
    labels = [f'{cov} bp, {n} intervals' 
              for cov, n in zip(df1_coverage_and_count.coverage, df1_coverage_and_count['count'])])

bf.vis.plot_intervals(df2, show_coords=True, xlim=(0,16), colors='lightpink')
```

### Complement
```{eval-rst}
Equally important to finding which genomic features overlap is finding those that do not. :func:`bioframe.complement` returns a BedFrame of intervals not covered by any intervals in an input BedFrame. 

Complments are defined relative to a `genomic view`, or ViewFrame. A genomic view is an ordered set of non-overlapping genomic intervals with a unique set of names, defined by a string ‘name’. See :ref:`Technical Notes<Technical_Notes>` for more details and :func:`bioframe.core.checks.is_viewframe` for implementation. `genomic views` can be provided in a number of formats, including a dict of {str:int} that defines chromosome name and length (see :func:`bioframe.core.construction.from_any` for details).
```

```{code-cell} ipython3
df_complemented = bf.complement(df1, view_df={'chr1':15})
display(df_complemented)
```

```{code-cell} ipython3
bf.vis.plot_intervals(df_complemented, show_coords=True, xlim=(0,16), colors='lightpink')
```

If no view is provided, complement is calculated per unique chromosome in the input with right limits of ``np.iinfo(np.int64).max``.
```{code-cell} ipython3
df_complemented = bf.complement(df1)
display(df_complemented)
```

### Flexible column naming

+++

Genomic analyses often deal with dataframes with inhomogeneously named columns. Bioframe offers a way to set the default column names that are most convenient for your analyses. 

Default bedframe column names are stored in ``bioframe.core.specs_rc``. 

```{code-cell} ipython3
bf.core.specs._rc
```

If the dataframes we wish to work with have `['CHROMOSOME', 'LEFT', 'RIGHT']`, we can either pass cols to operations in ``bioframe.ops``:

```{code-cell} ipython3
df1_diff_colnames = pd.DataFrame([
    ['chr1', 1, 5],
    ['chr1', 3, 8]],
    columns=['CHROMOSOME', 'LEFT', 'RIGHT']
)

df2_diff_colnames = pd.DataFrame([
    ['chr1', 4, 8],
    ['chr1', 10, 11]],
    columns=['CHROMOSOME', 'LEFT', 'RIGHT']
)
```

```{code-cell} ipython3
bf.overlap(
    df1_diff_colnames, df2_diff_colnames,
    cols1=['CHROMOSOME', 'LEFT', 'RIGHT'],
    cols2=['CHROMOSOME', 'LEFT', 'RIGHT'],
)
```

Or, we can update the default column names:

```{code-cell} ipython3
with bf.core.specs.update_default_colnames(['CHROMOSOME', 'LEFT', 'RIGHT']):
    display(bf.overlap(
        df1_diff_colnames, df2_diff_colnames,
    ))
```

```{code-cell} ipython3
# setting colnames back to default.
bf.core.specs.update_default_colnames(['chrom', 'start', 'end'])
bf.core.specs._rc
```

