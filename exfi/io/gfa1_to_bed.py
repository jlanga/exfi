#!/usr/bin/env python3

"""exfi.io.gfa1_to_bed.py: submodule to read a GFA1 file and convert it to BED4
"""

import logging

import numpy as np

from exfi.io.bed import BED4_COLS, BED4_DTYPES
from exfi.io.read_gfa import read_gfa1

def gfa1_to_bed4(filename):
    """Read a GFA1 file and convert it to BED4"""

    logging.info('Converting GFA1 to BED4')
    containments = read_gfa1(filename)['containments']

    logging.info('Renaming columns')
    containments = containments.rename(columns={
        "container": "chrom",
        "contained": "name"
    })
    logging.info('Overlap to int')
    containments["overlap"] = containments\
        .overlap.map(lambda x: np.int(x[:-1]))
    logging.info('Computing coordinates')
    containments["chrom_start"] = containments["pos"]
    containments["chrom_end"] = containments["pos"] + containments["overlap"]
    containments = containments[BED4_COLS]
    containments = containments.astype(BED4_DTYPES)
    logging.info('Done')
    return containments
