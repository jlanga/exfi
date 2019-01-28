#!/usr/bin/env python3

"""exfi.io.bed4_to_gfa1.py: submodule to write a BED4 dataframe to GFA1 format
"""

import logging

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
    logging.info('Computing the header')
    header = pd.DataFrame(
        data=[["H", "VN:Z:1.0"]],
        columns=HEADER_COLS
    )
    return header


def compute_segments(bed4, transcriptome_dict, masking='none'):
    """Create the Segments subdataframe for GFA1 file"""
    logging.info('Computing node2sequence')
    segments = bed4_to_node2sequence(
        bed4=bed4, transcriptome_dict=transcriptome_dict
    )
    logging.info('Computing edge2overlap')
    edge2overlap = bed4_to_edge2overlap(bed4)
    logging.info('Masking')
    segments = mask(
        node2sequence=segments, edge2overlap=edge2overlap, masking=masking
    )
    del edge2overlap

    # Add the S and length columns
    logging.info('Adding the record_type')
    segments["record_type"] = "S"

    return segments[SEGMENT_COLS]


def compute_links(bed4):
    """Compute the Links subdataframe of a GFA1 file."""
    logging.info('Computing edge2overlap')
    links = bed4_to_edge2overlap(bed4=bed4)\
        .rename(columns={'u': 'from', 'v': 'to'})
    logging.info('Adding record_type, from_orient, to_orient')
    links["record_type"] = "L"
    links["from_orient"] = "+"
    links["to_orient"] = "+"
    logging.info('Computing the overlap between exons')
    links["overlap"] = links.overlap.map(lambda x: str(x) + "M" if x >= 0 else str(-x) + "N")
    logging.info('Reordering')
    return links[LINK_COLS]


def compute_containments(bed4):
    """Create the minimal containments subdataframe"""
    containments = bed4.copy()
    logging.info('Adding record_type, container, container_orient, contained, '
                 'contained_orient, and pos')
    containments["record_type"] = "C"
    containments["container"] = containments["chrom"]
    containments["container_orient"] = "+"
    containments["contained"] = containments["name"]
    containments["contained_orient"] = "+"
    containments["pos"] = containments["chrom_start"]
    logging.info('Computing the overlap')
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
        logging.info('Writing the header')
        compute_header()\
            .to_csv(gfa, sep="\t", header=False, index=False)
    with open(gfa1_fn, "a", 1024**3) as gfa:
        logging.info('Writing the segments')
        compute_segments(
            bed4=bed4, transcriptome_dict=transcriptome_dict, masking=masking
            )\
            .to_csv(gfa, sep="\t", header=False, index=False)
        logging.info('Writing the links')
        compute_links(bed4=bed4)\
            .to_csv(gfa, sep="\t", header=False, index=False)
        logging.info('Writing the containments')
        compute_containments(bed4=bed4)\
            .to_csv(gfa, sep="\t", header=False, index=False)
        logging.info('Writing the paths')
        compute_paths(bed4=bed4)\
            .to_csv(gfa, sep="\t", header=False, index=False)
