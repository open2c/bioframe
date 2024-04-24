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

If you work with bed files you can simply load them using `read_table`, it will
create a pandas [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html)
which supports all the bioframe operations.

Altertantively if you want to work on `gtf` files, you do not need to turn them
into bed files, you can directly read them (with e.g. [gtfparse](https://github.com/openvax/gtfparse/tree/master))
and turn them into bedframe by renaming the `seqname` column into `chrom`.

Any DataFrame object with `'chrom'`, `'start'`, and `'end'` columns will support
all the following operations TODO `API_fileops`

## `bedtools intersect`

### Original unique entries from the first bed `-u`

Note that this gives one row per overlap and can contain duplicates,

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

### Original entries from the first bed for each overlap`-wa`

Note that this gives one row per overlap and can contain duplicates,

```sh
bedtools intersect -wa -a A.bed -b B.bed > out.bed
```

```py
overlap = bf.overlap(A, B, how='inner', suffixes=('_1','_2'), return_index=True)
out = A.loc[overlap['index_1']]
```

### Original entries from the second bed `-wb`

```sh
bedtools intersect -wb -a A.bed -b B.bed > out.bed
```

```py
overlap = bf.overlap(A, B, how='inner', suffixes=('_1','_2'), return_index=True)
out = B.loc[overlap['index_2']]
```

### Intersect with multiple beds

```sh
bedtools intersect -wa -a A.bed -b B.bed C.bed D.bed> out.bed
```

```py
others = pd.concat([B, C, D])
overlap = bf.overlap(A, others, how='inner', suffixes=('_1','_2'), return_index=True)
out = A.loc[overlap['index_1']]
```

### Keep no overlap `-v`

```sh
bedtools intersect -wa -a A.bed -b B.bed -v > out.bed
```

```py
out = bf.setdiff(A, B)
```

### Force strandedness `-s`

For intersection

```sh
bedtools intersect -wa -a A.bed -b B.bed -s > out.bed
```

```py
overlap = bf.overlap(A, B, on=['strand'], suffixes=('_1','_2'), return_index=True, how='inner')
out = A.loc[overlap['index_1']]
```

For non intersection

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
```

