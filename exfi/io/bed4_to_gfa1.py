#!/usr/bin/env python3

"""exfi.io.bed4_to_gfa1.py: submodule to write a BED4 dataframe to GFA1 format
"""

import pandas as pd

from exfi.io.bed import \
    bed4_to_node2sequence, \
    bed4_to_edge2overlap

from exfi.io.gfa1 import \
    HEADER_COLS, SEGMENT_COLS, LINK_COLS, CONTAINMENT_COLS, PATH_COLS

from exfi.io.masking import \
    mask

def compute_header():
    """Write GFA1 header"""
    header = pd.DataFrame(
        data=[["H", "VN:Z:1.0"]],
        columns=HEADER_COLS
    )
    return header


def compute_segments(bed4, transcriptome_dict, masking='none'):
    """Create the Segments subdataframe for GFA1 file"""
    segments = bed4_to_node2sequence(
        bed4=bed4, transcriptome_dict=transcriptome_dict
    )
    edge2overlap = bed4_to_edge2overlap(bed4)
    segments = mask(
        node2sequence=segments, edge2overlap=edge2overlap, masking=masking
    )
    del edge2overlap


    # Add the S and length columns
    segments["record_type"] = "S"

    # Compute lengths
    segments["length"] = segments\
        .sequence.map(lambda x: "LN:i:" + str(len(x)))

    # reorder
    segments = segments\
        [SEGMENT_COLS]

    return segments


def compute_links(bed4):
    """Compute the Links subdataframe of a GFA1 file."""
    links = bed4_to_edge2overlap(bed4=bed4)\
        .rename(columns={'u': 'from', 'v': 'to'})
    links["record_type"] = "L"
    links["from_orient"] = "+"
    links["to_orient"] = "+"
    links["overlap"] = links.overlap.map(lambda x: str(x) + "M" if x >= 0 else str(-x) + "N")
    links = links[LINK_COLS]
    return links


def compute_containments(bed4):
    """Create the minimal containments subdataframe"""
    containments = bed4.copy()
    containments["record_type"] = "C"
    containments["container"] = containments["chrom"]
    containments["container_orient"] = "+"
    containments["contained"] = containments["name"]
    containments["contained_orient"] = "+"
    containments["pos"] = containments["chrom_start"]
    containments["overlap"] = containments["chrom_end"] - containments["chrom_start"]
    containments["overlap"] = containments.overlap.map(lambda x: str(x) + "M")
    containments = containments.drop(
        ["chrom", "chrom_start", "chrom_end", "name"], axis=1
    )
    return containments[CONTAINMENT_COLS]


def compute_paths(bed4):
    """Compute the Paths section of the GFA1 file"""
    paths = bed4.copy()
    paths["name"] = paths["name"].map(lambda x: x + "+")
    paths = paths\
        .drop(columns=["chrom_start", "chrom_end"])\
        .groupby("chrom", axis=0)\
        .aggregate(lambda x: ",".join(x.tolist()))
    paths = paths.astype({"name": str})  # It may end up as float
    paths = paths.reset_index(drop=False)
    paths["record_type"] = "P"
    paths = paths.rename({"chrom": "path_name", "name": "segment_names"}, axis=1)
    paths["overlaps"] = "*"
    paths = paths[PATH_COLS]
    return paths


def bed4_to_gfa1(gfa1_fn, bed4, transcriptome_dict, masking='none'):
    """Convert the BED4 dataframe into a GFA1 file"""
    with open(gfa1_fn, "w", 1024**3) as gfa:
        compute_header()\
            .to_csv(gfa, sep="\t", header=False, index=False)
    with open(gfa1_fn, "a", 1024**3) as gfa:
        compute_segments(
            bed4=bed4, transcriptome_dict=transcriptome_dict, masking=masking
            )\
            .to_csv(gfa, sep="\t", header=False, index=False)
        compute_links(bed4=bed4)\
            .to_csv(gfa, sep="\t", header=False, index=False)
        compute_containments(bed4=bed4)\
            .to_csv(gfa, sep="\t", header=False, index=False)
        compute_paths(bed4=bed4)\
            .to_csv(gfa, sep="\t", header=False, index=False)
