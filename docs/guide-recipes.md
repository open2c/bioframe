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

## Obtain overlapping intervals with matching strandedness?
Use overlap with the ``on`` argument:
```
df = bf.overlap(df1, df2, on=[‘strand’])
```

## Obtain overlapping intervals with opposite strandedness?
Overlap then filter pairs of opposite strandedness:
```
df = bf.overlap(df1, df2)
df = df.loc[df["strand"]!=df["strand_"]]
```
## Obtain intervals that exceed 50% coverage by another set of intervals?
Coverage, then filter pairs by fractional coverage:
```
df = bf.coverage(df1, df2)
df = df[ ( df["coverage"] / (df["end"]-df["start"]) ) >=0.50]
```

## Shift all intervals on the positive strand by 10bp?
Use pandas indexing:
```
df.loc[df.strand=="+",["start", "end"]] += 10
```

## Obtain intervals overlapped by at least 2 intervals from another set?
Count overlaps, then filter:
```
df = bf.count_overlaps(df1, df2)
df = df[ df["count"] >= 2]
```

## Find strand-specific downstream genomic features?
Use closest after filtering by strand, and passing the `ignore_upsream=True` argument.
```
bioframe.closest(df1.loc[df1['strand']=='+'], df2, ignore_upstream=True)
```

For gener, the upstream/downstream direction might be defined by the direction of transcription.
Use `direction_col='strand'` to set up the direction:
```
bioframe.closest(df1, df2, ignore_upstream=True, direction_col='strand')
```

## Drop non-autosomes from a bedframe?
Use pandas DataFrame.isin(values):
```
df[ ~df.chrom.isin(['chrX','chrY'])]
```
