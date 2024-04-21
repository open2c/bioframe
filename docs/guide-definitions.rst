.. _Definitions:

Definitions
===========

Interval:
    - An *interval* is a tuple of integers (start, end) with start <= end.
    - Coordinates are assumed to be 0-based and intervals half-open (1-based ends) i.e. [start, end).
    - An interval has a *length* equal to (end - start).
    - A special case where start and end are the same, i.e. [X, X), is interpreted as a *point* (aka an *empty interval*, i.e. an edge between 1-bp bins). A point has zero length.
    - Negative coordinates are permissible for both ends of an interval.

Properties of a pair of intervals:
    - Two intervals can either *overlap*, or not. The overlap length = max(0, min(end1, end2) - max(start1, start2)). Empty intervals can have overlap length = 0.
    - When two intervals overlap, the shorter of the two intervals is said to be *contained* in the longer one if the length of their overlap equals the length of the shorter interval. This property is often referred to as nestedness, but we use the term “contained” as it is less ambiguous when describing the relationship of sets of intervals to one interval.
    - If two intervals do not overlap, they have a *distance* = max(0, max(start1, start2) - min(end1, end2)).
    - If two intervals have overlap=0 and distance=0, they are said to be *abutting*.

Scaffold:
    - A chromosome, contig or, more generally, a *scaffold* is an interval defined by a unique string and has a length>=0, with start=0 and end=length, implicitly defining an interval [0, length).

Genome assembly:
    - The complete set of scaffolds associated with a genome is called an *assembly* (e.g. defined by the reference sequence from NCBI, etc.).

Genomic interval:
    - A *genomic interval* is an interval with an associated scaffold, or chromosome, defined by a string, i.e. a triple (chrom, start, end).
    - Genomic intervals on different scaffolds never overlap and do not have a defined distance.
    - Genomic intervals can extend beyond their associated scaffold (e.g. with negative values or values greater than the scaffold length), as this can be useful in downstream applications. If they do, they are not contained by their associated scaffold.
    - A *base-pair* is a special case of a genomic interval with length=1, i.e. (chrom, start, start+1)
    - *strand* is an (optional) property of a genomic interval which specifies an interval’s orientation on its scaffold. Note start and end are still defined with respect to the scaffold’s reference orientation (positive strand), even if the interval lies on the negative strand. Intervals on different strands can either be allowed to overlap or not.

View (i.e. a set of Genomic Regions):
    - A genomic *view* is an ordered set of non-overlapping genomic intervals each having a unique name defined by a string. Individual named intervals in a view are *regions*, defined by a quadruple, e.g. (chrom, start, end, name).
    - A view thus specifies a unified 1D coordinate system, i.e. a projection of multiple genomic regions onto a single axis.
    - We define views separately from the scaffolds that make up a genome assembly, as a set of more constrained and ordered genomic regions are often useful for downstream analysis and visualization.
    - An assembly is a special case of a view, where the individual regions correspond to the assembly’s entire scaffolds.

Associating genomic intervals with views
    - Similarly to how genomic intervals are associated with a scaffold, they can also be associated with a region from a view with an additional string, making a quadruple (chrom, start, end, view_region). This string must be *cataloged* in the view, i.e. it must match the name of a region in the view. Typically the interval would be contained in its associated view region, or, at the minimum, have a greater overlap with that region than other view regions.
    - If each interval in a set is contained in their associated view region, the set is *contained* in the view.
    - A set of intervals *covers* a view if each region in the view is contained by the union of its associated intervals. Conversely, if a set does not cover all of view regions, the interval set will have *gaps* relative to that view (stretches of bases not covered by an interval).

Properties of sets of genomic intervals:
    - A set of genomic intervals may have overlaps or not. If it does not, it is said to be *overlap-free*.
    - A set of genomic intervals is *tiling* if it: (i) covers the associated view, (ii) is contained in that view, and (iii) is overlap-free. Equivalently, a tiling set of intervals (a) has an initial interval that begins at the start of each region and (b) a final interval that terminates at the end of each region, and (c) every base pair is associated with a unique interval.
