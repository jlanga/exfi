#!/usr/bin/env python3

"""exfi.io.read_bed.py: BED importer"""

import pandas as pd
import numpy as np


def read_bed3(filename):
    """Read a BED file and return the BED3 dataframe."""
    bed3 = pd.read_csv(
        filepath_or_buffer=filename,
        header=None,
        sep='\t',
        usecols=[0, 1, 2],
        names=["chrom", "chromStart", "chromEnd"],
        dtype={"chrom": np.str, "chromStart": np.int64, "chromEnd": np.int64},
        engine='c'
    )
    return bed3
