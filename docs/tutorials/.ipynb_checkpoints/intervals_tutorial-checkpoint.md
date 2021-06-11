---
jupytext:
  formats: ipynb,md:myst
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

Genomic interval operations 
===========================

```automodule:: bioframe.ops
---
autosummary:
members:
---
```

# Interval Operations

```{code-cell} ipython3
import itertools

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

import bioframe as bf
import bioframe.vis
```

### Example interval sets

```{code-cell} ipython3
df1 = pd.DataFrame([
    ['chr1', 1, 5],
    ['chr1', 3, 8],
    ['chr1', 8, 10],
    ['chr1', 12, 14]],
    columns=['chrom', 'start', 'end']
)

df2 = pd.DataFrame([
    ['chr1', 4, 8],
    ['chr1', 10, 11]],
    columns=['chrom', 'start', 'end']
)
```

```{code-cell} ipython3
bf.vis.plot_intervals(df1, show_coords=True, xlim=(0,16))
plt.title('set 1')

bf.vis.plot_intervals(df2, show_coords=True, xlim=(0,16), colors='lightpink')
plt.title('set 2')
```

### Overlap

```{code-cell} ipython3
overlapping_intervals = bf.overlap(df1, df2, how='inner')
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

```{code-cell} ipython3
df_annotated = bf.cluster(df1, min_dist=0)
display(df_annotated)
bf.vis.plot_intervals(df_annotated, labels=df_annotated['cluster'], show_coords=True, xlim=(0,16))
```

```{code-cell} ipython3
df_annotated = bf.cluster(df1, min_dist=None)
display(df_annotated)
bf.vis.plot_intervals(df_annotated, labels=df_annotated['cluster'], show_coords=True, xlim=(0,16))
```

```{code-cell} ipython3
df_annotated = bf.cluster(df1, min_dist=2)
display(df_annotated)
bf.vis.plot_intervals(df_annotated, labels=df_annotated['cluster'], show_coords=True, xlim=(0,16))
```

### Merge

```{code-cell} ipython3
df_merged = bf.merge(df1, min_dist=0)

display(df_merged)
bf.vis.plot_intervals(df_merged, show_coords=True, xlim=(0,16))
```

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
closest_intervals = bf.closest(df1, df2)
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
closest_intervals = bf.closest(df1, df2, ignore_overlaps=True)
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
n_overlaps = bf.overlap(df1, df2, return_index=True).groupby('index_1').agg({"index_2": "count"})
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

```{code-cell} ipython3

```
