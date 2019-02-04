#!/usr/bin/env python3

"""tests.compare.compare.py: Auxiliary variables for tests.compare"""

import pandas as pd

from exfi.io.read_bed import read_bed3
from exfi.io.bed import BED3_COLS, BED3_DTYPES
from exfi.compare import \
    TP_DF_COLS, TP_DF_DTYPES, \
    STATS_COLS, STATS_DTYPES


BED3_EMPTY_FN = "tests/compare/empty.bed"
BED3_TRUE_FN = "tests/compare/true.bed"
BED3_PRED_FN = "tests/compare/pred.bed"


BED3_EMPTY = read_bed3(BED3_EMPTY_FN)
BED3_TRUE = read_bed3(BED3_TRUE_FN)
BED3_PRED = read_bed3(BED3_PRED_FN)


TP_DF = pd.DataFrame(
    data=[
        ['ENSDART00000000004', '474', '542', 'ENSDART00000000004', '474', '542', '68'],
        ['ENSDART00000000004', '542', '1311', 'ENSDART00000000004', '542', '1311', '769'],
        ['ENSDART00000000004', '1413', '2475', 'ENSDART00000000004', '1413', '2476', '1062'],
        ['ENSDART00000000005', '0', '1655', 'ENSDART00000000005', '0', '1656', '1655'],
        ['ENSDART00000000005', '1656', '1812', 'ENSDART00000000005', '1656', '1813', '156'],
        ['ENSDART00000000005', '1812', '1949', 'ENSDART00000000005', '1813', '1950', '136'],
        ['ENSDART00000000005', '2289', '2603', 'ENSDART00000000005', '2290', '2604', '313']
    ]
)

FP_DF = pd.DataFrame(
    data=[
        ['ENSDART00000000004', '0', '47'],
        ['ENSDART00000000004', '45', '240'],
        ['ENSDART00000000004', '239', '339'],
        ['ENSDART00000000004', '338', '474'],
        ['ENSDART00000000004', '1310', '1414'],
        ['ENSDART00000000005', '1948', '2101'],
        ['ENSDART00000000005', '2100', '2207'],
        ['ENSDART00000000005', '2205', '2292']
    ]
)


FN_DF = pd.DataFrame(
    data=[
        ['ENSDART00000000004', '0', '48'],
        ['ENSDART00000000004', '48', '241'],
        ['ENSDART00000000004', '241', '340'],
        ['ENSDART00000000004', '340', '474'],
        ['ENSDART00000000004', '1311', '1413'],
        ['ENSDART00000000005', '1950', '2102'],
        ['ENSDART00000000005', '2102', '2206'],
        ['ENSDART00000000005', '2206', '2290']
    ]
)


TRUE_POSITIVES = TP_DF\
    .rename(columns={i: j for i, j in enumerate(TP_DF_COLS)})\
    .astype(dtype=TP_DF_DTYPES)\
    .drop(columns=6)

FALSE_POSITIVES = FP_DF\
    .rename(columns={i: j for i, j in enumerate(BED3_COLS)})\
    .astype(BED3_DTYPES)

FALSE_NEGATIVES = FN_DF\
    .rename(columns={i: j for i, j in enumerate(BED3_COLS)})\
    .astype(BED3_DTYPES)

CLASSIFICATION_EMPTY = {
    'true_positives': pd.DataFrame(columns=TP_DF_COLS)\
        .astype(dtype=TP_DF_DTYPES),
    'false_positives': pd.DataFrame(columns=BED3_COLS)\
        .astype(dtype=BED3_DTYPES),
    'false_negatives': pd.DataFrame(columns=BED3_COLS)\
        .astype(dtype=BED3_DTYPES)
}

CLASSIFICATION = {
    'true_positives': TRUE_POSITIVES,
    'false_positives': FALSE_POSITIVES,
    'false_negatives': FALSE_NEGATIVES
}

STATS_PER_EXON = pd.DataFrame(
    data=[[
        15.0, 14.0, 7.0, 8.0, 8.0,
        0.4666666666666667, 0.4666666666666667, 0.40651851851851856
    ]],
    columns=STATS_COLS
)\
.astype(STATS_DTYPES)


STATS_PER_BASE = pd.DataFrame(
    data=[[
        5080.0, 5090.0, 4159.0, 931.0, 911.0,
        0.8170923379174853, 0.8203155818540434, 2.1950225255006406
    ]],
    columns=STATS_COLS
)\
.astype(STATS_DTYPES)
