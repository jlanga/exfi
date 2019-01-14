#!/usr/bin/env python3

"""exfi.io.read_bed3.py: BED3 importer"""


def read_bed3(filename):
    """Read a BED file and return the BED3 dataframe."""
    import pandas as pd
    import numpy as np
    bed3 = pd.read_table(
        filepath_or_buffer=filename,
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "chromStart", "chromEnd"],
        dtype={"chrom": np.str, "chromStart": np.int, "chromEnd": np.int},
        engine='c'
    )
    return bed3
