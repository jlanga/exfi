#!/usr/bin/env python3

"""
Module to build the associated splice graph from a set of bed3 records
"""

import logging

import networkx as nx
import pandas as pd

def _bed3_to_str(bed3_record):
    """(str, int, int) -> str

    BED3 to string
    """
    if len(bed3_record) == 3:
        return "{0}:{1}-{2}".format(*bed3_record)
    else:
        raise IndexError("Incorrect number of elements in record")



def bed3_records_to_bed6df(iterable_of_bed3):
    """(iterable of tuples) -> pd.dataframe

    Convert an iterable of bed3 records to a bed6 as dataframe
    bed6 =[seqid, start, end, name, strand, score]
    """
    logging.info("\tbed3 -> bed6")
    bed6_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    return pd.DataFrame(
        data=(bed3_record + (_bed3_to_str((bed3_record)), 0, '+')
              for bed3_record in iterable_of_bed3),
        columns=bed6_cols
    )\
    .sort_values(by=bed6_cols[0:2])



def bed6df_to_node2coordinates(bed6df):
    """(pd.DataFrame) -> dict

    Get from the BED6 dataframe the correspondece of name -> (chrom, start, end)
    """
    logging.info("\tbed6df -> node2coordinates")
    # Check for extreme case:
    if bed6df.shape[0] > 0:

        # Compute the node_id: coordinates dict
        node2coordinate = bed6df\
            .sort_values(['chrom', 'start', 'end'])\
            .drop(["score", "strand"], axis=1)\
            .assign(
                coordinates=bed6df\
                    [["chrom", "start", "end"]]\
                    .apply(tuple, axis=1)
            )\
            .drop(["chrom", "start", "end"], axis=1)\
            .set_index("name", "coordinates")\
            .to_dict()["coordinates"]

        # Reprocess the dict, one node may be in multiple transcripts at once
        node2coordinate = {
            key: (value,)
            for key, value in node2coordinate.items()
        }

        return node2coordinate
    return {}



def bed6df_to_path2node(bed6df):
    """(pandas.df) -> {transcript_id: (node1, ..., nodeN)}

    Get a dict containing transcript_id to the tuple of node names that compose
    it in order, indicating the path."""
    logging.info("\tbed6df -> path2node")
    if bed6df.shape[0] > 0:
        return bed6df\
            .sort_values(['chrom', 'start', 'end'])\
            .drop(['start', 'end', 'strand', 'score'], axis=1)\
            .rename(columns={'chrom':'path'})\
            .groupby('path')\
            .agg(lambda x: tuple(x.tolist()))\
            .to_dict()["name"]
    return {}



def compute_edge_overlaps(splice_graph):
    """(nx.DiGraph) -> dict

    Get the overlap between connected exons:
    - Positive overlap means that they overlap that number of bases,
    - Zero that they occur next to each other
    - Negative that there is a gap in the transcriptome of that number of bases
    (one or multiple exons of length < kmer)
    Return dict {(str, str): int} (node1, node2, and overlap)
    Note: the splice graph must have already the nodes written with coordinates,
    and the edges alredy entered too.

    Hypothesis: node2coords.values should only hold one value
    """
    logging.info("\tComputing edge overlaps")

    #Init
    node2coords = nx.get_node_attributes(
        G=splice_graph,
        name='coordinates'
    )

    edge_overlaps = {edge: None for edge in splice_graph.edges()}

    for (node1, node2) in edge_overlaps.keys():
        node1_end = node2coords[node1][0][2]
        node2_start = node2coords[node2][0][1]

        # Overlap in bases, 0 means one next to the other, negative numbers a gap
        overlap = node1_end - node2_start
        edge_overlaps[(node1, node2)] = overlap

    return edge_overlaps



def build_splice_graph(positive_exons_bed):
    """(Iterable of (str, int, int)) -> Nx.DiGraph

    Build the splice_graph from a list of bed records

    splice_graph is a directed graph, whose nodes
        - are an identifier, the tuple in string format
        - whose attributes are
            - the cooridnates in (str, int, int) format
    and whose edges
        - are connected exons in any way
        - attributes are the overlap between them:
            - positive means there is an overlap of that number of bases
            - zero means no overlap
            - negative means a gap of that number of bases
    """
    logging.info("Building splice graph")
    # Initialize graph
    splice_graph = nx.DiGraph()

    # Process bed records
    bed6df = bed3_records_to_bed6df(positive_exons_bed)

    # Add nodes
    logging.info("Adding nodes")
    splice_graph.add_nodes_from(bed6df["name"].tolist())
    nx.set_node_attributes(  # Add coordinates
        G=splice_graph,
        name="coordinates",
        values=bed6df_to_node2coordinates(bed6df)
    )

    # Edges
    logging.info("Adding edges")
    transcript2path = bed6df_to_path2node(bed6df)
    for path in transcript2path.values():
        splice_graph.add_path(path)
    nx.set_edge_attributes(
        G=splice_graph,
        name='overlaps',
        values=compute_edge_overlaps(splice_graph)
    )

    return splice_graph
