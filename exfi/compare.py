#!/usr/bin/env python3

"""exfi.compare.py: submodule to compare BED, GFA and GFF3 files to measure exfi's performance"""

import sys
import logging
from subprocess import Popen, PIPE
from os import remove
from tempfile import mkstemp

import pandas as pd
import numpy as np

from exfi.io.bed import \
    BED3_COLS, BED3_DTYPES


TP_DF_COLS = [
    'chrom_pred', 'chrom_start_pred', 'chrom_end_pred',
    'chrom_true', 'chrom_start_true', 'chrom_end_true'
]
TP_DF_DTYPES = {
    'chrom_pred': object, 'chrom_start_pred': np.int, 'chrom_end_pred': np.int,
    'chrom_true': object, 'chrom_start_true': np.int, 'chrom_end_true': np.int
}

STATS_COLS = [
    'true', 'predicted', 'true_positives', 'false_positives', 'false_negatives',
    'precision', 'recall', 'f_1'
]

STATS_DTYPES = {
    'true': np.float, 'predicted': np.float, 'true_positives': np.float,
    'false_positives': np.float, 'false_negatives': np.float,
    'precision': np.float, 'recall': np.float, 'f_1': np.float
}

def bedtools_intersect(bed1_fn, bed2_fn, additional_flags=None):
    """Bedtools intersect wrapper

    bed1_fn: path to BED file
    bed2_fn: path to BED file
    additional_flags: list with additional flags to pass to bedtools in the shape of a list:
        ['-r', '-v', '-f', '0.95']

    returns a list of lists ()
    """
    if additional_flags is None:
        additional_flags = []
    command = ["bedtools", "intersect", "-a", bed1_fn, "-b", bed2_fn] + \
        additional_flags
    process = Popen(command, stdout=PIPE, stderr=PIPE)
    data = [
        line.decode().strip().split()
        for line in process.stdout.readlines()
    ]
    process.stdout.close()
    status_code = process.wait()
    if status_code != 0:
        sys.exit(
            "ERROR: something went wrong:\n" + \
            ",".join([x.decode() for x in process.stderr.readlines()])
        )
    return data



def classify(bed3_true, bed3_pred, fraction=0.95):
    """Compute the True Posivites, False Positives and False Negatives
    between two BED3 dataframes (true and predicted) using bedtools intersect

    Fraction is a real number between 0 and 1 that stands for
    minimum similarity.

    Result is a dict where values are dataframes
    """
    logging.info('Classifying')

    bed3_pred_fn = mkstemp()[1]
    bed3_true_fn = mkstemp()[1]

    # Dump to disk
    logging.info("Dumping predictions to disk")
    bed3_true\
        .sort_values(by=['chrom', 'chrom_start', 'chrom_end'])\
        .to_csv(path_or_buf=bed3_true_fn, sep='\t', index=False, header=False)

    logging.info('Dumping true values to disk')
    bed3_pred\
        .sort_values(by=['chrom', 'chrom_start', 'chrom_end'])\
        .to_csv(path_or_buf=bed3_pred_fn, sep='\t', index=False, header=False)

    logging.info('Computing true positives')
    true_positives_df = pd.DataFrame(
        data=bedtools_intersect(
            bed1_fn=bed3_pred_fn,
            bed2_fn=bed3_true_fn,
            additional_flags=['-f', f'{fraction}', '-r', '-wo']
        ),
        columns=TP_DF_COLS + [6]
    ).astype(dtype=TP_DF_DTYPES)\
    .drop(columns=6)

    logging.info('Computing false positives')
    false_positives_df = pd.DataFrame(
        data=bedtools_intersect(
            bed1_fn=bed3_pred_fn,
            bed2_fn=bed3_true_fn,
            additional_flags=['-f', f'{fraction}', '-r', '-v'],
        ),
        columns=BED3_COLS
    ).astype(BED3_DTYPES)

    logging.info('Computing false negatives')
    false_negatives_df = pd.DataFrame(
        data=bedtools_intersect(
            bed1_fn=bed3_true_fn,
            bed2_fn=bed3_pred_fn,
            additional_flags=['-f', f'{fraction}', '-r', '-v']
        ),
        columns=BED3_COLS
    ).astype(BED3_DTYPES)

    remove(bed3_pred_fn)
    remove(bed3_true_fn)

    return {
        'true_positives': true_positives_df,
        'false_positives': false_positives_df,
        'false_negatives': false_negatives_df
    }



def compute_precision(true_positives, false_positives):
    """Compute precision

    >>> compute_precision(0, 10)
    0.0
    >>> compute_precision(446579, 13932)
    0.969747
    """
    return true_positives / (true_positives + false_positives)



def compute_recall(true_positives, false_negatives):
    """Compute recall

    >>> compute_recall(0, 10)
    0.0
    >>> compute_recall(446579, 48621)
    0.901815
    """
    return true_positives / (true_positives + false_negatives)



def compute_f_1(true_positives, false_positives, false_negatives):
    """Compute F_1

    >>> compute_f_1(0, 10, 10)
    0
    >>> compute_f_1(446579, 13932, 48621)
    0.934548
    """
    precision = compute_precision(true_positives, false_positives)
    recall = compute_recall(true_positives, false_negatives)
    return 2 * precision * recall / (precision + recall)



def compute_stats_per_exon(classification):
    """Compute the classification stats per exon

    Input should be the one from exfi.compare.classify
    """

    logging.info('Computing the stats per exon')

    tp_exons = classification['true_positives'].shape[0]
    fp_exons = classification['false_positives'].shape[0]
    fn_exons = classification['false_negatives'].shape[0]

    true_exons = tp_exons + fn_exons
    pred_exons = tp_exons + fp_exons

    stats = pd.DataFrame(
        data=[[
            true_exons, pred_exons,
            tp_exons, fp_exons, fn_exons,
            compute_precision(tp_exons, fp_exons),
            compute_recall(tp_exons, fn_exons),
            compute_f_1(tp_exons, fp_exons, fn_exons)
        ]],
        columns=STATS_COLS
    )\
    .astype(STATS_DTYPES)

    return stats



def compute_true_bases(classification):
    """Compute the total number of bases in the truth splice graph"""
    true_positives, _, false_negatives = classification.values()
    return \
        np.sum(true_positives.chrom_end_true - true_positives.chrom_start_true) + \
        np.sum(false_negatives.chrom_end - false_negatives.chrom_start)



def compute_pred_bases(classification):
    """Compute the total number of bases in the predicted splice graph"""
    true_positives, false_positives, _ = classification.values()
    return \
        np.sum(true_positives.chrom_end_pred - true_positives.chrom_start_pred) + \
        np.sum(false_positives.chrom_end - false_positives.chrom_start)



def compute_true_positive_bases(classification):
    """Compute the total number of true positive bases in the predicted splice
    graph

    TP bases are the common part between the two exons
    """
    true_positives = classification['true_positives']
    return np.sum(
        true_positives[['chrom_end_true', 'chrom_end_pred']].min(axis=1) -
        true_positives[['chrom_start_true', 'chrom_start_pred']].max(axis=1)
    )



def compute_false_positive_bases(classification):
    """Compute the total number of true positive bases in the predicted splice
    graph

    FP bases are:
        - all bases in false positive exons
        - all bases over predicted in the start
        - all bases over predicted in the end
    """
    tp_df, fp_df, _ = classification.values()

    starts = tp_df[tp_df.chrom_start_pred < tp_df.chrom_start_true]
    ends = tp_df[tp_df.chrom_end_pred > tp_df.chrom_end_true]

    fp_bases = \
        np.sum(fp_df.chrom_end - fp_df.chrom_start) + \
        np.sum(starts.chrom_start_true - starts.chrom_start_pred) + \
        np.sum(ends.chrom_end_pred - ends.chrom_end_true)

    return fp_bases



def compute_false_negative_bases(classification):
    """Compute the total number of true positive bases in the predicted splice
    graph

    FN bases are:
        - all bases in false negative exons
        - all bases underpredicted in the start
        - all bases underpredicted in the end
    """
    tp_df, _, fn_df = classification.values()

    starts = tp_df[tp_df.chrom_start_pred > tp_df.chrom_start_true]
    ends = tp_df[tp_df.chrom_end_pred < tp_df.chrom_end_true]

    fn_bases = \
        np.sum(fn_df.chrom_end - fn_df.chrom_start) + \
        np.sum(starts.chrom_start_pred - starts.chrom_start_true) + \
        np.sum(ends.chrom_end_pred - ends.chrom_end_true)

    return fn_bases



def compute_stats_per_base(classification):
    """Compute the classification stats per base pair

    Input should be the one from exfi.compare.classify
    """
    logging.info('Computing the stats per base')

    true_bases = compute_true_bases(classification)
    pred_bases = compute_pred_bases(classification)

    tp_bases = compute_true_positive_bases(classification)
    fp_bases = compute_false_positive_bases(classification)
    fn_bases = compute_false_negative_bases(classification)


    stats = pd.DataFrame(
        data=[[
            true_bases, pred_bases,
            tp_bases, fp_bases, fn_bases,
            compute_precision(tp_bases, fp_bases),
            compute_recall(tp_bases, fn_bases),
            compute_f_1(tp_bases, fp_bases, fn_bases)
        ]],
        columns=STATS_COLS
    )\
    .astype(STATS_DTYPES)

    return stats
