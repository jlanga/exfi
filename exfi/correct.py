#!/usr/bin/env python3

"""exfi.new_correct.py: fill small overlaps and gaps with abyss-sealer"""

import logging

from tempfile import \
    mkstemp

from subprocess import Popen
import os

import pandas as pd

from exfi.io.bed import \
    BED3_COLS, \
    bed3_to_bed4, \
    bed4_to_node2sequence, \
    bed4_to_edge2overlap


def prepare_sealer(bed4, transcriptome_dict, args):
    """exfi.new_correct.prepare_sealer: inspect the bed4 file and create a fasta
    file where pairs of exons have a small gap between them or have a small
    overlap.
    """

    sealer_input = mkstemp()

    max_fp_bases = args["max_fp_bases"]
    max_gap_size = args["max_gap_size"]

    node2sequence = bed4_to_node2sequence(bed4, transcriptome_dict)
    edge2overlap = bed4_to_edge2overlap(bed4)
    node2sequence_dict = node2sequence.set_index("name").to_dict()["sequence"]

    # Disable warnings
    pd.options.mode.chained_assignment = None

    # Compute the small gaps
    small_gaps = edge2overlap\
        .loc[(edge2overlap.overlap < 0) & (edge2overlap.overlap <= max_gap_size)]

    small_gaps["identifier"] = small_gaps['u'] + "~" + small_gaps['v']

    small_gaps["data_to_map"] = tuple(zip(small_gaps.u, small_gaps.v))

    small_gaps["sequence"] = small_gaps.data_to_map\
        .map(
            lambda x: \
            node2sequence_dict[x[0]][0:-max_fp_bases] + \
            100 * 'N' + \
            node2sequence_dict[x[1]][max_fp_bases:]
        )

    small_gaps = small_gaps[["identifier", "sequence"]]

    # Compute pairs of overlapping exons
    overlaps = edge2overlap.loc[edge2overlap.overlap >= 0]
    overlaps["data_to_map"] = tuple(zip(overlaps.u, overlaps.v, overlaps.overlap))
    overlaps["identifier"] = overlaps.u + "~" + overlaps.v
    overlaps["sequence"] = overlaps.data_to_map\
        .map(
            lambda x: \
            node2sequence_dict[x[0]][0:-x[2] - max_fp_bases] + \
            100 * 'N' + \
            node2sequence_dict[x[1]][x[2] + max_fp_bases:]
        )
    overlaps = overlaps[["identifier", "sequence"]]

    # Put again the warning
    pd.options.mode.chained_assignment = 'warn'

    # Merge the results
    for_sealer = pd.concat([small_gaps, overlaps])
    for_sealer["fasta"] = ">" + for_sealer["identifier"] + "\n" + for_sealer["sequence"] + "\n"
    for_sealer = for_sealer[["fasta"]]

    with open(sealer_input[1], "w", 1*1024**3) as f_in:
        for fasta_record in for_sealer.fasta.values:
            f_in.write(fasta_record)

    return sealer_input[1]


def run_sealer(sealer_input_fn: str, args: dict) -> str:
    """Run abyss-sealer with the parameters in args, and the scaffold in
    sealer_input.

    args = {
        "kmer": int,
        "max_gap_size": int,
        "input_bloom": str,
    }

    :param str sealer_input_fn: Input filename for sealer (the scaffold).
    :param dict args: Dict of argumnets for sealer

    """
    #logging.debug("\tRunning abyss-sealer")
    # Run sealer
    sealer_output_prefix = mkstemp()
    c_sealer = [
        'abyss-sealer',
        '--input-scaffold', sealer_input_fn,
        '--flank-length', str(args["kmer"]),
        '--max-gap-length', "30",
        '--kmer', str(args["kmer"]),
        '--fix-errors',
        '--input-bloom', args["bloom"],
        '--mask',
        '--output-prefix', sealer_output_prefix[1],
        '--verbose'
    ]

    # Execute
    p_sealer = Popen(c_sealer)
    p_sealer.communicate()

    # Clean files
    os.remove(sealer_output_prefix[1] + "_log.txt")
    os.remove(sealer_output_prefix[1] + "_scaffold.fa")
    os.remove(sealer_output_prefix[1])

    return sealer_output_prefix[1] + "_merged.fa"


def collect_sealer_results(filename):
    """Read the fasta output from sealer and return the merged nodes"""

    if os.path.getsize(filename) == 0:
        return pd.DataFrame(data=None, columns=["u", "v"])

    headers = pd.read_csv(filename, header=None, sep="\t")
    headers = headers.iloc[::2]  # Take odd rows: headers.

    headers.columns = ["raw"]
    headers["clean"] = headers\
        .raw\
        .str.slice(1)\
        .str.rsplit("_", 2).str[0]\
        .str.split("~")
    headers["u"], headers["v"] = headers.clean.str
    headers = headers[["u", "v"]]
    headers = headers.reset_index(drop=True)
    return headers


def apply_correction_to_bed4(bed4, sealed_edges):
    """Merge nodes into a single ones, being careful with the coordinates"""
    if sealed_edges.shape[0] == 0:
        return bed4
    new_bed4 = bed4.copy().set_index("name")
    for row in sealed_edges.iloc[::-1].itertuples():
        new_bed4.loc[row.u, "chrom_end"] = new_bed4.loc[row.v, "chrom_end"]
    new_bed4 = new_bed4.drop(sealed_edges.v)
    new_bed4 = bed3_to_bed4(new_bed4[BED3_COLS])
    return new_bed4.reset_index(drop=True)


def correct_bed4(bed4, transcriptome_dict, args):
    """Inspect the bed4 for small gaps and overlaps, write a fasta file for
    sealer, and correct the bed4.
    """
    logging.info('Preparing abyss-sealer')
    sealer_input_fn = prepare_sealer(
        bed4=bed4, transcriptome_dict=transcriptome_dict, args=args
    )
    logging.info('Running abyss-sealer')
    sealer_output_fn = run_sealer(sealer_input_fn=sealer_input_fn, args=args)
    logging.info('Collecting abyss-sealer\'s results')
    sealer_results = collect_sealer_results(filename=sealer_output_fn)
    logging.info('Applying correction to BED4')
    bed4_corrected = apply_correction_to_bed4(bed4, sealer_results)
    os.remove(sealer_input_fn)
    os.remove(sealer_output_fn)
    return bed4_corrected
