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

# How do I

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
overlap = bf.overlap(A, B, how='inner', suffixes=('_1','_2'), return_index=True)
out = A.loc[~A.index.isin(set(overlap['index_1'].unique()))]
```

