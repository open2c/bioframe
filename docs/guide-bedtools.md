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

# Emulating bedtools commands

If you want to work on `gtf` files, you do not need to turn them into bed files,
you can directly read them (with e.g. [gtfparse](https://github.com/openvax/gtfparse/tree/master))
and turn them into bedframe by renaming the `seqname` column into `chrom`.

## `bedtools intersect`

### Original entries from the first bed `-wa`

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
overlap = bf.overlap(A, B, on=['strand'], suffixes=('_1','_2'), return_index=True)
out = A.loc[overlap['index_1']]
```

For non intersection

```sh
bedtools intersect -wa -a A.bed -b B.bed -v -s > out.bed
```

```py
out = bf.setdiff(A, B, on=['strand'])
```

### Minimum overlap a as fraction of A `-f`

We want to keep rows of A that are covered at least 70% by elements from B

```sh
bedtools intersect -wa -a A.bed -b B.bed -f 0.7 > out.bed
```

```py
cov = bf.coverage(A, B)
out = A[cov['coverage'] / (cov['end']-cov['start']) ) >=0.70]
```
