Technical Notes
===============

        
Interval:
    - A tuple of integers (start,end) with start <= end.  
    - Coordinates are assumed to be 0-based and intervals half-open (1-based ends) i.e. [start, end). A special case of [X, X) is interpreted as a point (empty interval, i.e. an edge between 1-bp bins). 
    - An interval has a length equal to (end - start). An empty interval has zero length. 
    - Negative coordinates are permissible for both ends.
    
Properties of pairs of intervals:
    - can either have an overlap, or not. Overlap = max(0, min(end1, end2) - max(start1, start2)). Empty intervals can have overlap = 0.
    - if they overlap, are they contained. They are contained if the length of their overlap equals the length of the shorter of the two intervals. This property is often referred to as nested, but we use contained as it is less ambiguous for sets of intervals relative to one interval. 
    - if they do not overlap, they have a distance = max(0, max(start1, start2) - min(end1, end2)). 
    
Scaffold:
    - a chromosome, contig or, more generally, a scaffold is defined by a unique string and has a length>=0, with start=0 and end=length, implicitly defining an interval [0, length).
    - a scaffold defines an independent coordinate axis upon which intervals can be defined.
    
Genome: 
    - the complete set of scaffolds associated with a genome is called an assembly (e.g. defined by the reference sequence from NCBI, etc.).
    
Genomic interval:
    - an interval with an associated scaffold, or chromosome, defined by a string ‘chrom’. A genomic interval is thus a triple (chrom, start, end). 
    - strand is an (optional) binary property of a genomic interval which specifies an interval’s orientation on its scaffold. Note start and end are still defined with respect to the scaffold’s reference orientation (positive strand), even if the interval lies on the negative strand.
    - intervals on different scaffolds never overlap and do not have a defined distance. Intervals on different strands can either be allowed to overlap or not.
    - note that genomic intervals can extend beyond their associated scaffold (e.g. with negative values or values greater than the scaffold length, for use in downstream applications).
    - a base-pair is a special case of a genomic interval with length=1, i.e. (chrom, start, start+1)
    
View (i.e. a set of Genomic Regions):
    - a genomic view is an ordered set of non-overlapping genomic intervals with a unique set of names, defined by a string ‘name’. Individual members of this set are regions and are quadruples, (chrom, start, end, name). A view thus specifies a global coordinate axis, i.e. a conversion from a genomic coordinate system to a single axis.
    - we define a view separately from scaffolds, as a set of more constrained genomic intervals are often useful for downstream analysis, including: pileups, P(s), saddles, dotfinding, and analyses by chromosome arm.
    - any genomic interval can be associated with a view via an additional string. This string must be cataloged in the view, i.e. it must match the name of a region in the view, and is by default defined by ‘view_region.’ Note if no ‘view_region’ is assigned to a genomic interval, then ‘chrom’ implicitly defines an associated region
    - an assembly is a special case of a set of genomic regions, where scaffolds are individual regions.

Properties of sets of genomic intervals:
    - if all intervals are contained in their associated view region, the set is contained.
    - if the set of all associated regions are contained in the union of the set of intervals, the set of intervals covers the set of regions. Conversely, if a set does not cover all of its regions, the interval set will have gaps (stretches of bases not covered by an interval).
    - may have overlaps or not. If it does not, it is said to be overlap-free. An overlap-free interval set can be depicted as a stepwise-constant function.
    - a set of intervals is tiling if it: (i) covers all regions in the set of associated regions, (ii) is contained, and (iii) is overlap-free. Equivalently, a tiling set of intervals (a) has an initial interval that begins at the start of each region and (b) a final interval that terminates at the end of each region, and (c) every base pair is associated with a unique interval. 
    - the coverage depth at each base pair in a view is defined as the number of genomic intervals overlapping that base pair. Thus, when an interval set possesses overlaps, one or more bases will have an interval coverage depth greater than 1. Conversely, a tiling set has a uniform coverage depth of 1. Similarly, gaps have a coverage depth of 0.
    - Any overlap-free, contained set of genomic intervals (e.g. a bedGraph) can thus be converted to a tiling set of genomic intervals by filling the gaps with explicit intervals.

Properties of BedFrames (i.e. dataframes of genomic intervals):
    - three required columns determine the sets of genomic intervals (default = (chrom, start, end))
    - other reserved but not required column names: (strand, name, view_region).

        - name is expected to be unique
        - view_region is expected to point to a region in a view with a matching name
        - strand is expected to be encoded with strings (‘+’, ‘-’, ‘.’). 

    - additional columns are allowed: zodiac_sign, soundcloud, twitter_name, etc.
    - repeated intervals are allowed
    - meant to be compatible with BED format
    - pandas index is not used
    - intervals can be sorted:

        - scaffold-sorted: on (chrom, start, end), where chr is sorted lexicographically
        - view-sorted: on (view_region, start, end) where view_region is sorted by order in the view
        
    - null values are allowed, but only as pd.NA (using np.nan is discouraged as it results in unwanted type re-casting)

ViewFrame:
    - genomic interval dataframe where:

        - items  are a quadruple (chrom, start, end, name)
        - items in name column are unique (each region has a unique name), and by default provides a lookup table for the view_region
        - no null values
        - intervals are non-overlapping

    - views can be defined by: 
        
        - dictionary of string:ints (start=0 assumed) or string:tuples (start,end)
        - pandas series chromsizes (start=0, name=chrom)

    - often a (sorted) set of contigs



