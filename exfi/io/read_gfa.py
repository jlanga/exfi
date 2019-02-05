#!/usr/bin/env python3

"""exfi.io.read_gfa.py: submodule to read GFA1 files"""

import logging
import pandas as pd

from exfi.io.gfa1 import \
    HEADER_COLS, SEGMENT_COLS, LINK_COLS, CONTAINMENT_COLS, PATH_COLS, \
    HEADER_DTYPES, SEGMENT_DTYPES, LINK_DTYPES, CONTAINMENT_DTYPES, PATH_DTYPES


def _process_header(data):
    """Process the header from raw data"""
    logging.info('Processing header lines')
    return pd.DataFrame(
        data=[x[0:2] for x in data if x[0] == 'H'],
        columns=HEADER_COLS
    ).astype(HEADER_DTYPES)


def _process_segments(data):
    """Process the segments from raw data"""
    logging.info('Processing segment lines')
    return pd.DataFrame(
        data=[x[0:3] for x in data if x[0] == "S"],
        columns=SEGMENT_COLS
    ).astype(SEGMENT_DTYPES)


def _process_links(data):
    """"Process the links from raw data"""
    logging.info('Processing link lines')
    return pd.DataFrame(
        data=[x[0:6] for x in data if x[0] == 'L'],
        columns=LINK_COLS
    ).astype(LINK_DTYPES)


def _process_containments(data):
    """Process the containments from raw data"""
    logging.info('Processing containment lines')
    return pd.DataFrame(
        data=[x[0:7] for x in data if x[0] == 'C'],
        columns=CONTAINMENT_COLS
    ).astype(CONTAINMENT_DTYPES)


def _process_paths(data):
    """Process the paths from raw data"""
    logging.info('Processing path lines')
    return pd.DataFrame(
        data=[x[0:4] for x in data if x[0] == 'P'],
        columns=PATH_COLS
    ).astype(PATH_DTYPES)


def read_gfa1(gfa1_fn):
    """Read the GFA1 file in gfa1_fn and return a dict of dataframes where the
    keys are header, segments, links, containments, and paths. Values are
    DataFrames, with the exception of the header"""

    logging.info('Reading a GFA1 from disk')

    with open(gfa1_fn, 'r') as gfa:

        gfa1 = {}

        logging.info('Reading raw data')
        data = [
            x.strip().split("\t")
            for x in gfa.readlines() if x[0] in set(['H', 'S', 'L', 'C', 'P'])
        ]

        gfa1['header'] = _process_header(data)
        gfa1['segments'] = _process_segments(data)
        gfa1['links'] = _process_links(data)
        gfa1['containments'] = _process_containments(data)
        gfa1['paths'] = _process_paths(data)


        logging.info('Done')
        return gfa1
