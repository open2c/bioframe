import numpy as np
import pandas as pd
import pytest

import bioframe


def test_to_ucsc_colorstring():
    bioframe.to_ucsc_colorstring("red") == "255,0,0"
    bioframe.to_ucsc_colorstring("blue") == "0,0,255"
    bioframe.to_ucsc_colorstring("green") == "0,128,0"
    bioframe.to_ucsc_colorstring("black") == "0,0,0"
    bioframe.to_ucsc_colorstring("white") == "255,255,255"
    bioframe.to_ucsc_colorstring("r") == "255,0,0"
    bioframe.to_ucsc_colorstring("tomato") == "255,99,71"
    bioframe.to_ucsc_colorstring("xkcd:sky blue") == "135,206,235"
    bioframe.to_ucsc_colorstring("#abc") == "170,187,204"
    bioframe.to_ucsc_colorstring("#ff0000") == "255,0,0"
    bioframe.to_ucsc_colorstring("#ff000055") == "255,0,0"
    bioframe.to_ucsc_colorstring((1, 0, 0)) == "255,0,0"
    bioframe.to_ucsc_colorstring((1, 0, 0, 0.5)) == "255,0,0"
    bioframe.to_ucsc_colorstring((0, 0, 1)) == "0,0,255"
    bioframe.to_ucsc_colorstring(None) == "0"
    bioframe.to_ucsc_colorstring("none") == "0"
    bioframe.to_ucsc_colorstring(np.nan) == "0"
    bioframe.to_ucsc_colorstring(pd.NA) == "0"

    with pytest.raises(ValueError):
        bioframe.to_ucsc_colorstring("notacolor")

    df = bioframe.from_any(
        [
            ["chr1", 0, 10, "red"],
            ["chr1", 10, 20, "blue"],
            ["chr2", 0, 10, "green"],
            ["chr2", 10, 20, None],
        ]
    )
    df["itemRgb"] = df["name"].apply(bioframe.to_ucsc_colorstring)
    assert df["itemRgb"].tolist() == ["255,0,0", "0,0,255", "0,128,0", "0"]
