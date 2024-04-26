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

# Bioframe for bedtools users


Bioframe is built around the analysis of genomic intervals as a pandas [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html) in memory, rather than working with tab-delimited text files saved on disk.

Bioframe supports reading a number of standard genomics text file formats via [`read_table`](https://bioframe.readthedocs.io/en/latest/api-fileops.html#bioframe.io.fileops.read_table), including BED files (see [schemas](https://github.com/open2c/bioframe/blob/main/bioframe/io/schemas.py)), which will load them as pandas DataFrames, a complete list of helper functions is [available here](API_fileops).

Any DataFrame object with `'chrom'`, `'start'`, and `'end'` columns will support the genomic [interval operations in bioframe](API_ops). The names of these columns can also be customized via the `cols=` arguments in bioframe functions.

For example, with gtf files, you do not need to turn them into bed files, you can directly read them into pandas (with e.g. [gtfparse](https://github.com/openvax/gtfparse/tree/master)). For gtfs, it is often convenient to rename the `'seqname'` column to `'chrom'`, the default column name used in bioframe.

Finally, if needed, bioframe provides a convenience function to write dataframes to a standard BED file using [`to_bed`](https://bioframe.readthedocs.io/en/latest/api-fileops.html#bioframe.io.bed.to_bed).


## `bedtools intersect`

### Select unique entries from the first bed overlapping the second bed `-u`

```sh
bedtools intersect -u -a A.bed -b B.bed > out.bed
```

```py
overlap = bf.overlap(A, B, how='inner', suffixes=('_1','_2'), return_index=True)
out = A.loc[overlap['index_1'].unique()]
```

### Report the number of hits in B `-c`

Reports 0 for A entries that have no overlap with B.

```sh
bedtools intersect -c -a A.bed -b B.bed > out.bed
```

```py
out = bf.count_overlaps(A, B)
```

### Return entries from both beds for each overlap `-wa -wb`

```sh
bedtools intersect -wa -wb -a A.bed -b B.bed > out.bed
```

```py
out = bf.overlap(A, B, how='inner')
```

**Note:** This is called an "inner join", and is analogous to an inner pandas join or merge. The default column suffixes in the output dataframe are `''` (nothing) for A's columns and `'_'` for B's columns.

### Include all entries from the first bed, even if no overlap `-loj`

```sh
bedtools intersect -wa -wb -loj -a A.bed -b B.bed > out.bed
```

```py
out = bf.overlap(A, B, how='left')
```

**Note:** This is called a "left-outer join".

### Select entries from the first bed for each overlap `-wa`

```sh
bedtools intersect -wa -a A.bed -b B.bed > out.bed
```

```py
overlap = bf.overlap(A, B, how='inner', suffixes=('_1','_2'), return_index=True)
out = A.loc[overlap['index_1']]

# Alternatively
out = bf.overlap(A, B, how='inner')[A.columns]
```

> **Note:** This gives one row per overlap and can contain duplicates. The output dataframe of the former method will use the same pandas index as the input dataframe `A`, while the latter result --- the join output --- will have an integer range index, like a pandas merge.

### Select entries from the second bed for each overlap `-wb`

```sh
bedtools intersect -wb -a A.bed -b B.bed > out.bed
```

```py
overlap = bf.overlap(A, B, how='inner', suffixes=('_1','_2'), return_index=True)
out = B.loc[overlap['index_2']]

# Alternatively
out = bf.overlap(A, B, how='inner', suffixes=('_', ''))[B.columns]
```

> **Note:** This gives one row per overlap and can contain duplicates. The output dataframe of the former method will use the same pandas index as the input dataframe `B`, while the latter result --- the join output --- will have an integer range index, like a pandas merge.


### Intersect multiple beds against A

```sh
bedtools intersect -wa -a A.bed -b B.bed C.bed D.bed > out.bed
```

```py
others = pd.concat([B, C, D])
overlap = bf.overlap(A, others, how='inner', suffixes=('_1','_2'), return_index=True)
out = A.loc[overlap['index_1']]
```

### Return everything in A that doesn't overlap with B `-v`

```sh
bedtools intersect -wa -a A.bed -b B.bed -v > out.bed
```

```py
out = bf.setdiff(A, B)
```

**Note:** We call this a set difference.

### Force strandedness `-s`

For intersection

```sh
bedtools intersect -wa -a A.bed -b B.bed -s > out.bed
```

```py
overlap = bf.overlap(A, B, on=['strand'], suffixes=('_1','_2'), return_index=True, how='inner')
out = A.loc[overlap['index_1']]
```

For non-intersection `-v`

```sh
bedtools intersect -wa -a A.bed -b B.bed -v -s > out.bed
```

```py
out = bf.setdiff(A, B, on=['strand'])
```

### Minimum overlap as a fraction of A `-f`

We want to keep rows of A that are covered at least 70% by elements from B

```sh
bedtools intersect -wa -a A.bed -b B.bed -f 0.7 > out.bed
```

```py
cov = bf.coverage(A, B)
out = A.loc[cov['coverage'] / (cov['end'] - cov['start']) ) >= 0.70]

# Alternatively
out = bf.coverage(A, B).query('coverage / (end - start) >= 0.7')[A.columns]
```
