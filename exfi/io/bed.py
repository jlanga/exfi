#!/usr/bin/env python3

"""exfi.io.bed.py: submodule to wrangle BED dataframes"""

import numpy as np

BED3_COLS = ['chrom', 'chrom_start', 'chrom_end']
BED3_DTYPES = {'chrom': np.str, 'chrom_start': np.int64, 'chrom_end': np.int64}


BED4_COLS = BED3_COLS + ['name']
BED4_DTYPES = {
    'chrom': np.str,
    'chrom_start': np.int64,
    'chrom_end': np.int64,
    'name': np.str
}



def bed3_to_bed4(bed3):
    """Take a BED3 dataframe and add the name as:
        "chrom:chromStart+chromEnd"
    """
    bed4 = bed3.copy()
    bed4["name"] = \
        bed4.chrom + ":" + \
        bed4.chrom_start.map(str) + "-" + \
        bed4.chrom_end.map(str)
    return bed4.sort_values(['chrom', 'chrom_start'])


def bed4_to_node2coordinates(bed4):
    """Compute the node2coordinates DataFrame: exon name, chrom, start, end"""
    node2coordinates = bed4\
        [["name", "chrom", "chrom_start", "chrom_end"]]\
        .set_index("name")
    return node2coordinates


def bed4_to_path2nodes(bed4):
    """Compute the correspondance between a transcript and its exons:
    {transcript_id : list of exons}.
    """
    return bed4\
        .drop(columns=["chrom_start", "chrom_end"])\
        .groupby("chrom")\
        .agg(lambda x: x.tolist())\
        .to_dict()["name"]


def bed4_to_node2sequence(bed4, transcriptome_dict):
    """Compute the correspondence between an exon name and its sequence as a
    pd.DataFrame with cols name and sequence.
    """
    node2sequence = bed4.copy()
    node2sequence["sequence"] = node2sequence.chrom.map(transcriptome_dict)
    node2sequence["data_to_map"] = list(zip(
        node2sequence.sequence,
        node2sequence.chrom_start,
        node2sequence.chrom_end
    ))
    node2sequence.sequence = node2sequence.data_to_map.map(lambda x: x[0][x[1]:x[2]])
    return node2sequence[["name", "sequence"]]


def bed4_to_edge2overlap(bed4):
    """Compute the overlaps between a pair of overlapping exons.
    Dataframe is name_u, name_v, overlap_int.
    """
    overlaps = bed4.copy()
    # Get the transcript_id of the next exon
    overlaps["chrom_next"] = overlaps["chrom"].shift(-1)
    # Get the name of the next exon
    overlaps["name_next"] = overlaps["name"].shift(-1)
    # Get the start of the next exon
    overlaps["chrom_start_next"] = overlaps["chrom_start"].shift(-1)
    # Get the end of the next exon
    overlaps["chrom_end_next"] = overlaps["chrom_end"].shift(-1)
    # Remove rows with different transcripts
    overlaps = overlaps\
        [overlaps["chrom"] == overlaps["chrom_next"]]
    # Convert types
    overlaps = overlaps.astype({"chrom_start_next": int, "chrom_end_next": int})
    # Compute the overlap
    overlaps["overlap"] = overlaps["chrom_end"] - overlaps["chrom_start_next"]
    # Convert again just in case
    overlaps.astype({"overlap": int})
    # Select and rename
    overlaps = overlaps\
        [["name", "name_next", "overlap"]]\
        .rename({"name": "u", "name_next": "v"}, axis=1)
    return overlaps.reset_index(drop=True)
