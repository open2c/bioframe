# Release notes

## [v0.3.1](https://github.com/open2c/bioframe/compare/v0.3.0...v0.3.1)

Date : 2021-11-15

API changes:

* `bioframe.sort_bedframe` does not append columns or modify their dtypes.


## [v0.3.0](https://github.com/open2c/bioframe/compare/v0.2.0...v0.3.0)

Date : 2021-08-31

Conceptual changes:
* we formulated strict definitions for genomic intervals, dataframes, and 
    their various properties. All bioframe functions are expected to follow
    to these definitions tightly.  

API changes:
* reorganize modules: 
    * ops - operations on genomic interval dataframes 
    * extras - miscellaneous operations, most involving
        genomic sequences and gene annotations
    * vis - visualizations of genomic interval dataframes
    * core.arrops - operations on genomic interval arrays
    * core.checks - tests for definitions of genomic interval dataframes
    * core.construction - construction and sanitation of genomic interval dataframes
    * core.specs - specifications for the implementation of genomic intervals in pandas.dataframes 
        (i.e. column names, datatypes, etc)
    * core.stringops - operations on genomic interval strings
    * io.fileops - I/O on common file formats for genomic data
    * io.schemas - schemas for standard tabular formats for genomic data storage
    * io.resources - interfaces to popular online genomic data resources 

* new functions: extras.pair_by_distance, ops.sort_bedframe, ops.assign_view, 
    dataframe constructors

* existing functions:
    * expand: take negative values and fractional values
    * overlap: change default suffixes, keep_order=True
    * subtract: add return_index and keep_order

* enable pd.NA for missing values, typecasting

New data:
* add schemas for bedpe, gap, UCSCmRNA, pgsnp
* add tables with curated detailed genome assembly information

Bugfixes:
* None?..

Miscellaneous:
* speed up frac_gc is faster now
* drop support for Python 3.6, add support for 3.9


## [v0.2.0](https://github.com/open2c/bioframe/compare/v0.1.0...v0.2.0)

Date : 2020-12-02

API changes
* `read_chromsizes` and `fetch_chromsizes`: add new `as_bed` parameter.
* `read_chromsizes` and `fetch_chromsizes`: revert to filtering chromosome names by default, but clearly expose `filter_chroms` kwarg.

Bug fixes
* Fixed `bioframe.split`
* Restored `frac_genome_coverage`


## [v0.1.0](https://github.com/open2c/bioframe/compare/v0.0.12...v0.1.0)

Date : 2020-09-23

First beta release.

### What's new

* New extensive dataframe genomic interval arithmetic toolsuite.
* Improved region handling and region querying functions.
* [Documentation!](https://bioframe.readthedocs.io/)

### Maintenance

* Dropped Python 2 support
* Refactoring of various genome operations and resources.
* Improved testing and linting
