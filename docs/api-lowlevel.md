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

# Low-level API

```{eval-rst}
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   lowlevel/arrops
   lowlevel/specs
   lowlevel/stringops

```

Low level array-based operations are used to implement the genomic interval operations on dataframes.

```{code-cell} ipython3
import itertools

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

import bioframe as bf
import bioframe.vis


from bioframe.core import arrops
```

```{code-cell} ipython3
starts1, ends1 = np.array([
    [1,5],
    [3,8],
    [8,10],
    [12,14]
]).T

starts2, ends2 = np.array([
    [4,8],
    [10,11],
]).T
```

```{code-cell} ipython3
bf.vis.plot_intervals_arr(
    starts = starts1,
    ends = ends1,
    xlim = (-0.5,14.5),
    labels = np.arange(0,starts1.shape[0]),
    show_coords = True)

bf.vis.plot_intervals_arr(
    starts = starts2,
    ends = ends2,
    colors = 'lightpink',
    xlim = (-0.5,14.5),
    labels = np.arange(0,starts2.shape[0]),
    show_coords = True)
```

```{code-cell} ipython3
arrops.overlap_intervals(starts1, ends1, starts2, ends2)
```

```{code-cell} ipython3
arrops.overlap_intervals_outer(starts1, ends1, starts2, ends2)
```

```{code-cell} ipython3
arrops.merge_intervals(starts1, ends1, min_dist=0)
```

```{code-cell} ipython3
arrops.merge_intervals(starts1, ends1, min_dist=None)
```

```{code-cell} ipython3
arrops.merge_intervals(starts1, ends1, min_dist=2)
```

```{code-cell} ipython3
arrops.complement_intervals(starts1, ends1)
```
