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
## Obtain intervals that exceed 50% overlap?
Overlap, then filter pairs on percentage overlap:
```
df = bf.overlap(df1, df2, return_overlap=True)
df = df[ (df["overlap_end"]-df["overlap_start"]) / (df["end"]-df["start"]) >=0.50]
```

## Shift all intervals on the positive strand by 10bp?
Use pandas indexing: 
```
df.loc[df.strand=="+",["start", "end"]] += 10
```




