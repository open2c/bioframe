.. _Specifications:

Specifications
===========

BedFrame (i.e. genomic intervals stored in a pandas dataframe):
    - In a BedFrame, three required columns specify the set of genomic intervals (default column names = (‘chrom’, ‘start’, ‘end’)).
    - Other reserved but not required column names: (‘strand’, ‘name’, ‘view_region’).

        - entries in column ‘name’ are expected to be unique
        - ‘view_region’ is expected to point to an associated region in a view with a matching name
        - ‘strand’ is expected to be encoded with strings (‘+’, ‘-’, ‘.’).

    - Additional columns are allowed: ‘zodiac_sign’, ‘soundcloud’, ‘twitter_name’, etc.
    - Repeated intervals are allowed.
    - The native pandas DataFrame index is not intended to be used as an immutable lookup table for genomic intervals in BedFrame. This is because many common genomic interval operations change the number of intervals stored in a BedFrame.
   - Two useful sorting schemes for BedFrames are:

        - scaffold-sorted: on (chrom, start, end), where chrom is sorted lexicographically.
        - view-sorted: on (view_region, start, end) where view_region is sorted by order in the view.

    - Null values are allowed, but only as pd.NA (using np.nan is discouraged as it results in unwanted type re-casting).
   - Note if no ‘view_region’ is assigned to a genomic interval, then ‘chrom’ implicitly defines an associated region
    - Note the BedFrame specification is a natural extension of the BED format ( ​​https://samtools.github.io/hts-specs/BEDv1.pdf ) for pandas DataFrames.

ViewFrames (a genomic view stored in a pandas dataframe)
    - BedFrame where:

           - intervals are non-overlapping
           - “name” column is mandatory and contains a set of unique strings.

    - Note that a ViewFrame can potentially be indexed by the name column to serve as a lookup table. This functionality is currently not implemented, because within the current Pandas implementation indexing by a column removes the column from the table.
    - Note that views can be defined by:

        - dictionary of string:ints (start=0 assumed) or string:tuples (start,end)
        - pandas series of chromsizes (start=0, name=chrom)
