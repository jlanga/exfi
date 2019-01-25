#!/usr/bin/env python3

"""exfi.io.read_gfa.py: submodule to read GFA1 files"""


import pandas as pd

from exfi.io.gfa1 import \
    HEADER_COLS, SEGMENT_COLS, LINK_COLS, CONTAINMENT_COLS, PATH_COLS, \
    HEADER_DTYPES, SEGMENT_DTYPES, LINK_DTYPES, CONTAINMENT_DTYPES, PATH_DTYPES


def read_gfa1(gfa1_fn):
    """Read the GFA1 file in gfa1_fn and return a dict of dataframes where the
    keys are header, segments, links, containments, and paths. Values are
    DataFrames, with the exception of the header"""

    with open(gfa1_fn, 'r') as gfa:

        gfa1 = {}

        data = [
            x.strip().split("\t")
            for x in gfa.readlines() if x[0] in set(['H', 'S', 'L', 'C', 'P'])
        ]

        gfa1['header'] = pd.DataFrame(
            data=[x[0:2] for x in data if x[0] == 'H'],
            columns=HEADER_COLS
        ).astype(HEADER_DTYPES)

        gfa1['segments'] = pd.DataFrame(
            data=[x[0:3] for x in data if x[0] == "S"],
            columns=SEGMENT_COLS
        ).astype(SEGMENT_DTYPES)

        gfa1['links'] = pd.DataFrame(
            data=[x[0:6] for x in data if x[0] == 'L'],
            columns=LINK_COLS
        ).astype(LINK_DTYPES)

        gfa1['containments'] = pd.DataFrame(
            data=[x[0:7] for x in data if x[0] == 'C'],
            columns=CONTAINMENT_COLS
        ).astype(CONTAINMENT_DTYPES)

        gfa1['paths'] = pd.DataFrame(
            data=[x[0:4] for x in data if x[0] == 'P'],
            columns=PATH_COLS
        ).astype(PATH_DTYPES)

        return gfa1
