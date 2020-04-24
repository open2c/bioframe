import pandas as pd
import bioframe


def test_expand():
    fake_bioframe = pd.DataFrame(
        {"chrom": ["chr1", "chr1", "chr2"], "start": [1, 50, 100], "end": [5, 55, 200]}
    )
    fake_chromsizes = {"chr1": 60, "chr2": 300}
    expand_bp = 10
    fake_expanded = bioframe.expand(fake_bioframe.copy(), expand_bp, fake_chromsizes)
    print(fake_expanded)
    assert fake_expanded.iloc[0].start == 0  # don't expand below zero
    assert (
        fake_expanded.iloc[1].end == fake_chromsizes["chr1"]
    )  # don't expand above chromsize
    assert (
        fake_expanded.iloc[2].end == fake_bioframe.iloc[2].end + expand_bp
    )  # expand end normally
    assert (
        fake_expanded.iloc[2].start == fake_bioframe.iloc[2].start - expand_bp
    )  # expand start normally
