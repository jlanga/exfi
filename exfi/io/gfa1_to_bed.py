#!/usr/bin/env python3

"""exfi.io.gfa1_to_bed.py: submodule to read a GFA1 file and convert it to BED4
"""

import numpy as np

from exfi.io.bed import BED4_COLS, BED4_DTYPES
from exfi.io.read_gfa import read_gfa1

def gfa1_to_bed4(filename):
    """Read a GFA1 file and convert it to BED4"""

    containments = read_gfa1(filename)['containments']

    containments = containments.rename(columns={
        "container": "chrom",
        "contained": "name"
    })
    containments["overlap"] = containments\
        .overlap.map(lambda x: np.int(x[:-1]))
    containments["chrom_start"] = containments["pos"]
    containments["chrom_end"] = containments["pos"] + containments["overlap"]
    containments = containments[BED4_COLS]
    containments = containments.astype(BED4_DTYPES)
    return containments
