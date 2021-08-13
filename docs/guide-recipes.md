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






