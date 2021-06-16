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

### DataFrames of genomic intervals

```{eval-rst}
The core object in bioframe are genomic interval dataframes, or bedframes. These can either be defined with :py:class:`pandas.DataFrame`:
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

BedFrames satisfy the following properties:  
- chrom, start, end columns  
- columns have valid dtypes (object/string/categorical, int, int)  
- all starts < ends.  
```{eval-rst}
Whether a dataframe satisfies these properties can be checked with :func:`bioframe.core.checks.is_bedframe`:
```
```{code-cell} ipython3
bioframe.is_bedframe(df2)
```
```{eval-rst}
See the :ref:`Technical Notes<Technical_Notes>` for more details.
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

### Overlaps
The two dataframes defined above contain two pairs of overlapping intervals:
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

### Cluster
```{eval-rst}
To demonstrate the usage of :py:func:`bioframe.cluster`, we use the same df1 as above:
```
bioframe.cluster
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

Extending the minimum distance to 2 makes all intervals part of the same cluster "0":
```{code-cell} ipython3
df_annotated = bf.cluster(df1, min_dist=2)
display(df_annotated)
bf.vis.plot_intervals(df_annotated, labels=df_annotated['cluster'], show_coords=True, xlim=(0,16))
```

### Merge


```{eval-rst}
Instead of returning cluster assignments, :py:func:`bioframe.merge` returns a new dataframe of merged genomic intervals. 

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

### Complement

```{code-cell} ipython3
bf.vis.plot_intervals(df1, show_coords=True, xlim=(0,16))
```

```{code-cell} ipython3
df_complemented = bf.complement(df1)
display(df_complemented)
```

```{code-cell} ipython3
df_complemented = bf.complement(df1, view_df={'chr1':16})
display(df_complemented)
bf.vis.plot_intervals(df_complemented, show_coords=True, xlim=(0,16), colors='lightpink')
```

### Closest

```{code-cell} ipython3
bf.closest(df1, df2)
```

```{code-cell} ipython3
bf.closest(df1, df2, return_input=2)
```

```{code-cell} ipython3
closest_intervals = bf.closest(df1, df2, suffixes=('_1','_2'))
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

```{code-cell} ipython3
bf.closest(df1, None, k=2)
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

```{code-cell} ipython3
bf.closest(df1, None)
```

### Coverage

```{code-cell} ipython3
df1_coverage = bf.coverage(df1, df2)
df1_coverage
```

### Overlaps

```{code-cell} ipython3
n_overlaps = bf.overlap(df1, df2, return_index=True, suffixes=('_1','_2')).groupby('index_1').agg({"index_2": "count"})
df1_n_overlaps = ( pd.Series(np.zeros(df1.shape[0], dtype=np.int64), index=df1.index)
                    .add(n_overlaps["index_2"], fill_value=0)
                    .astype(np.int64) )
df1_coverage['n_overlaps'] = df1_n_overlaps
df1_coverage
```

```{code-cell} ipython3
bf.vis.plot_intervals(
    df1_coverage, 
    show_coords=True, xlim=(0,16), 
    labels = [f'{cov} bp, {n} intervals' 
              for cov, n in zip(df1_coverage.coverage, df1_coverage['n_overlaps'])])

bf.vis.plot_intervals(df2, show_coords=True, xlim=(0,16), colors='lightpink')
```

### Column names

+++

Genomic analyses often deal with dataframes with inhomogeneously named columns. Bioframe offers a way to set the default column names that are most convenient for your analyses. 

Default bedframe column names are stored in bioframe.core.specs_rc. 

```{code-cell} ipython3
bf.core.specs._rc
```

If the dataframes we wish to work with have `['CHROMOSOME', 'LEFT', 'RIGHT']`, we can either pass cols to operations in bioframe.ops:

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

