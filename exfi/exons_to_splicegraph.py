#!/usr/bin/env python3

from exfi.exons_to_gapped_transcript import \
    build_transcript_to_exon_dict
from Bio import SeqIO
import networkx as nx
import pandas as pd


def exons_to_df(exons):
    """Convert an indexed fasta (exons) into a dataframe (~BED6)"""
    results = []
    for exon in exons.values():
        exon_id = exon.id
        transcript_coords = exon.description.split(" ")[1:]
        for transcript_coord in transcript_coords:
            transcript_id, coords = transcript_coord.split(":")
            start, end = coords.split("-")
            start = int(start)
            end = int(end)
            results.append(
                (transcript_id, start, end, exon_id, 0, "+")
            )
    return pd.DataFrame(
            data=results,
            columns=['transcript_id', 'start', 'end', 'exon_id', 'score', 'strand']
        )\
        .sort_values(
            by=['transcript_id', 'start','end']
        )


def exon_to_coordinates(exons_index):
    """Convert an indexed fasta (SeqIO.index) into a dict {exon_id : (transcript_id, start,
    end)} (str, int, int)"""
    exon_to_coord = {}
    for exon in exons_index.values():
        exon_id = exon.id
        # Drop id from desc
        transcript_coords = exon.description.split(" ")[1:]
        for transcript_coord in transcript_coords:
            # Compute values
            transcript_id, coords = transcript_coord.split(":")
            start, end = coords.split("-")
            start = int(start)
            end = int(end)
            # Add data
            if exon_id not in exon_to_coord:
                exon_to_coord[exon_id] = []
            exon_to_coord[exon_id].append((transcript_id, start, end))
    return exon_to_coord


def transcript_to_path(exon_df):
    """Get a Df containing transcript_id to list of exons, indicating the path"""
    return exon_df\
        .sort_values(['transcript_id', 'start', 'end'])\
        .drop(['start','end','score','strand'], axis=1)\
        .groupby('transcript_id')\
        .agg(lambda exon: exon.tolist())\
        .rename(columns={'exon_id':'path'})


def compute_edge_overlaps(splice_graph):
    """Get the overlap between connected exons:
    - Positive overlap means that they overlap that number of bases,
    - Zero that they occur next to each other
    - Negative that there is a gap in the transcriptome of that number of bases (one or multiple exons of length < kmer)

    Note: the splice graph must have already the nodes written with coordinates, and the edges alredy entered too.
    """
    #Init
    edge_overlaps = {}
    exon2coord = nx.get_node_attributes(
        G=splice_graph,
        name='coordinates'
    )

    for edge in splice_graph.edges():

        # Get involved nodes
        node1, node2 = edge

        # Get the list of transcripts that they belong
        node1_transcripts = set(coordinate[0] for coordinate  in exon2coord[node1])
        node2_transcripts = set(coordinate[0] for coordinate  in exon2coord[node2])
        intersection = node1_transcripts & node2_transcripts
        a_common_transcript = intersection.pop()

        # Get the end the first
        node1_coords = exon2coord[node1]
        node1_coords_in_transcript = [x for x in node1_coords if x[0] == a_common_transcript][0]
        node1_end = node1_coords_in_transcript[2]

        # Get the start of the next
        node2_coords = exon2coord[node2]
        node2_coords_in_transcript = [x for x in node2_coords if x[0] == a_common_transcript][0]
        node2_start = node2_coords_in_transcript[1]

        # Overlap in bases, 0 means one next to the other, negative numbers a gap
        overlap = node1_end - node2_start
        edge_overlaps[edge] = overlap

    return edge_overlaps