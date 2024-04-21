Resources
=========

Genome assembly metadata
------------------------

Bioframe provides a collection of genome assembly metadata for commonly used
genomes. These are accessible through a convenient dataclass interface via :func:`bioframe.assembly_info`.

The assemblies are listed in a manifest YAML file, and each assembly
has a mandatory companion file called `seqinfo` that contains the sequence
names, lengths, and other information. The records in the manifest file contain
the following fields:

- ``organism``: the organism name
- ``provider``: the genome assembly provider (e.g, ucsc, ncbi)
- ``provider_build``: the genome assembly build name (e.g., hg19, GRCh37)
- ``release_year``: the year of the assembly release
- ``seqinfo``: path to the seqinfo file
- ``cytobands``: path to the cytoband file, if available
- ``default_roles``: default molecular roles to include from the seqinfo file
- ``default_units``: default assembly units to include from the seqinfo file
- ``url``: URL to where the corresponding sequence files can be downloaded

The `seqinfo` file is a TSV file with the following columns (with header):

- ``name``: canonical sequence name
- ``length``: sequence length
- ``role``: role of the sequence or scaffold (e.g., "assembled", "unlocalized", "unplaced")
- ``molecule``: name of the molecule that the sequence belongs to, if placed
- ``unit``: assembly unit of the chromosome (e.g., "primary", "non-nuclear", "decoy")
- ``aliases``: comma-separated list of aliases for the sequence name

We currently do not include sequences with "alt" or "patch" roles in `seqinfo` files, but we
do support the inclusion of additional decoy sequences (as used by so-called NGS *analysis
sets* for human genome assemblies) by marking them as members of a "decoy" assembly unit.

The `cytoband` file is an optional TSV file with the following columns (with header):

- ``chrom``: chromosome name
- ``start``: start position
- ``end``: end position
- ``band``: cytogenetic coordinate (name of the band)
- ``stain``: Giesma stain result

The order of the sequences in the `seqinfo` file is treated as canonical.
The ordering of the chromosomes in the `cytobands` file should match the order
of the chromosomes in the `seqinfo` file.

The manifest and companion files are stored in the ``bioframe/io/data`` directory.
New assemblies can be requested by opening an issue on GitHub or by submitting a pull request.

.. automodule:: bioframe.io.assembly
   :autosummary:
   :members:

.. autoclass:: bioframe.io.assembly.GenomeAssembly
   :members:
   :undoc-members:


Remote resources
----------------
These functions now default to using the local data store, but can be used to obtain chromsizes or
centromere positions from UCSC by setting ``provider="ucsc"``.

.. automodule:: bioframe.io.resources
   :autosummary:
   :members:
