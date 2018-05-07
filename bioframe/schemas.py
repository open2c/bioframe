#-*- coding: utf-8 -*-

# https://genome.ucsc.edu/FAQ/FAQformat.html
BED12_FIELDS = ['chrom', 'start', 'end',
                'name', 'score', 'strand',
                'thickStart', 'thickEnd', 'rgb',
                'blockCount', 'blockSizes', 'blockStarts']
BED_FIELDS = BED12_FIELDS[:6]

BEDGRAPH_FIELDS = ['chrom', 'start', 'end', 'value']

BEDPE_FIELDS = ['chrom1', 'start1', 'end1',
                'chrom2', 'start2', 'end2',
                'name', 'score', 'strand1', 'strand2']
GFF_FIELDS = ['chrom', 'source', 'feature', 'start', 'end',
              'score', 'strand', 'frame', 'attributes']

PGSNP_FIELDS = ['chrom', 'start', 'end', 'name',
                'alleleCount', 'alleleFreq', 'alleleScores']
BEDRNAELEMENTS_FIELDS = ['chrom', 'start', 'end', 'name', 'score', 'strand',
              'level', 'signif', 'score2']
NARROWPEAK_FIELDS = ['chrom', 'start', 'end', 'name', 'score', 'strand',
                     'fc', '-log10p', '-log10q', 'relSummit']
BROADPEAK_FIELDS = ['chrom', 'start', 'end', 'name', 'score', 'strand',
                     'fc', '-log10p', '-log10q']
GAPPEDPEAK_FIELDS = ['chrom', 'start', 'end', 'name', 'score', 'strand',
                     'thickStart', 'thickEnd', 'rgb',
                     'blockCount', 'blockSizes', 'blockStarts',
                     'fc', '-log10p', '-log10q']
GAP_FIELDS = ['bin', 
              'chrom', 'start', 'end', 
              'ix', 'n', 
              'length', 'type', 'bridge']

UCSC_MRNA_FIELDS = [
    'bin', 'matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 
    'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand', 'qName', 'qSize', 
    'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 
    'blockSizes', 'qStarts', 'tStarts' ]
# http://ga4gh.org/#/fileformats-team
BAM_FIELDS = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR',
              'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'TAGs']
VCF_FIELDS = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

SCHEMAS = {
    'bed': BED_FIELDS,
    'bed3': BED_FIELDS[:3],
    'bed4': BED_FIELDS[:4],
    'bedGraph': BEDGRAPH_FIELDS,
    'bed5': BED_FIELDS[:5],
    'bed6': BED_FIELDS,
    'bed9': BED12_FIELDS[:9],
    'bed12': BED12_FIELDS,
    'gff': GFF_FIELDS,
    'gtf': GFF_FIELDS,
    'bedRnaElements': BEDRNAELEMENTS_FIELDS,
    'narrowPeak': NARROWPEAK_FIELDS,
    'broadPeak': BROADPEAK_FIELDS,
    'gappedPeak': GAPPEDPEAK_FIELDS,
    'bam': BAM_FIELDS,
    'vcf': VCF_FIELDS,
}

CHROM_NAME_PATTERNS = {
    'hg19': (r'^chr[0-9]+$', r'^chr[XY]$', r'^chrM$', r'^chr(?!Un).+_.*_random$', r'^chrUn_.*$', r'^chr(?!Un).+_.*_hap\d+$'),
    'hg38': (r'^chr[0-9]+$', r'^chr[XY]$', r'^chrM$', r'^chr(?!Un).+_.*_random$', r'^chrUn_.*$', r'^chr(?!Un).+_.*_alt$'),
    'mm9': (r'^chr[0-9]+$', r'^chr[XY]$', r'^chrM$', r'^chr(?!Un).+_random$', r'^chrUn_random$'),
    'mm10': (r'^chr[0-9]+$', r'^chr[XY]$', r'^chrM$', r'^chr(?!Un).+_random$', r'^chrUn_.*$'),
    'galGal4': (r'^chr[0-9]+$', r'^chr[ZW]$', r'^chrM$', r'^chrLGE64|chrLGE22C19W28_E50C23$',  r'^chr(?!Un).+_random$', r'^chrUn_.*$'),
    'galGal5': (r'^chr[0-9]+$', r'^chr[ZW]$', r'^chrM$', r'^chrLGE64$', r'^chr(?!Un).+_random$', r'^chrUn_.*$'),
    'dm3': (r'^chr[234][LR[Het]*]*$', r'^chr[XY][Het]*$', r'^chrM$', r'^chr[U][extra]*$'),
    'dm6': (r'^chr[234][LR]*', r'^chr[XY]$', r'^chrM$', r'^chr(?!Un).+_random$', r'^chrUn_.*$'),
    'ce10': (r'chr[IV]+$', r'^chrX$', r'^chrM$'),
    'ce11': (r'chr[IV]+$', r'^chrX$', r'^chrM$'),
    'sacCer3': (r'chr[IXV]+$', r'^chrM$'),
    'TAIR10': (r'^\d+', r'^MT|Pltd$')
}


UCSC_AUTOSOMES = r'^chr[0-9]+[A-Za-z]*$'
UCSC_SEXCHROMS = r'^chr[XYZW][A-Za-z]*$'
UCSC_NONNUCLEAR = r'^chrM$'
UCSC_OTHER = r'^chrLGE.*$'
UCSC_UNLOCALIZED = r'^chr(?!Un).+.*_random$'
UCSC_UNPLACED = r'^chrUn_.*$|^chrU[A-Za-z]*$'
UCSC_ALTCHROMS = r'^chr(?!Un).+_.*_hap\d+$|^chr(?!Un).+_.*_alt$'
NCBI_AUTOSOMES = r'^[0-9]+$'
NCBI_SEXCHROMS = r'^[XYZW]$'
NCBI_NONNUCLEAR = r'^MT$|^Pltd$'
ROMAN_LT10 = [r'^chrI+$', r'^chrIV$', r'^chrVI*$', r'^chrIX$']
ROMAN_LT20 = ROMAN_LT10 + [r'^chrX$', r'^chrXI*$', r'^chrXIV$', r'^chrXVI*$', r'^chrXIX$']
