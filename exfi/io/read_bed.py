#!/usr/bin/env python3

"""exfi.io.read_bed.py: BED importer"""

import pandas as pd

from exfi.io.bed import BED3_COLS, BED3_DTYPES


def read_bed3(filename):
    """Read a BED file and return the BED3 dataframe."""
    bed3 = pd.read_csv(
        filepath_or_buffer=filename,
        header=None,
        sep='\t',
        usecols=[0, 1, 2],
        names=BED3_COLS,
        engine='c'
    ).astype(BED3_DTYPES)
    return bed3
