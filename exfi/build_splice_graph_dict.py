#!/usr/bin/env python3

"""Module to build the associated splice graph from a set of bed3 records"""

from typing import \
    Iterable

import logging

import networkx as nx
import pandas as pd
import pathos.multiprocessing as mp

from exfi.classes import Node2Coordinates, Edge2Overlap, Coordinate, SpliceGraph, SpliceGraphDict, \
    Path2Nodes

def _bed3_to_str(bed3_record: Coordinate) -> str:
    """Convert a three element tuple into a string of the form a[0]:a[1]-a[2]

    :param tuple bed3_record: tuple. Must have length 3
    """
    if len(bed3_record) == 3:
        return "{0}:{1}-{2}".format(*bed3_record)
    else:
        raise IndexError("Incorrect number of elements in record")


def bed3_records_to_bed6df_dict(iterable_of_bed3: Iterable[Coordinate]) -> pd.DataFrame:
    """Convert an iterable of bed3 records to a bed6 as dataframe

    Columns of the bed6df are seqid, start, end, name, strand and score

    :param iterable_of_bed3: Iterable of tuples in BED3 format (seqid, start, end).
    """
    logging.info("\tbed3_records_to_bed6df_dict")
    bed6_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    bed6_df = pd.DataFrame(
        data=(bed3_record + (_bed3_to_str((bed3_record)), 0, '+')
              for bed3_record in iterable_of_bed3),
        columns=bed6_cols
    )\
        .sort_values(by=bed6_cols[0:2])\

    bed6df_dict = {
        transcript_id: dataframe
        for transcript_id, dataframe in bed6_df.groupby("chrom")
    }

    return bed6df_dict


def bed6df_to_node2coordinates(bed6df: pd.DataFrame) -> Node2Coordinates:
    """Get from the BED6 dataframe the dict of node_id : coordinates

    :param pd.DataFrame bed6df: DataFrame of BED6 records.
    """
    logging.debug("\tbed6bed6df_to_node2coordinates")
    # Check for extreme case:
    if bed6df.shape[0] == 0:
        return Node2Coordinates()
    # Compute the node_id: coordinates dict
    node2coordinate = bed6df\
        .sort_values(['chrom', 'start', 'end'])\
        .drop(["score", "strand"], axis=1)\
        .assign(
            coordinates=bed6df
            [["chrom", "start", "end"]]
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


def bed6df_to_path2node(bed6df: pd.DataFrame) -> Path2Nodes:
    """Get a dict containing transcript_id to the tuple of node names that compose it in order,
    indicating the path.

    :param pd.DataFrame bed6df: DataFrame with BED6 records.
    """
    logging.debug("\tbed6df -> path2node")
    if bed6df.shape[0] > 0:
        return bed6df\
            .sort_values(['chrom', 'start', 'end'])\
            .drop(['start', 'end', 'strand', 'score'], axis=1)\
            .rename(columns={'chrom': 'path'})\
            .groupby('path')\
            .agg(lambda x: tuple(x.tolist()))\
            .to_dict()["name"]
    return Path2Nodes()


def compute_edge_overlaps(splice_graph: SpliceGraph) -> Edge2Overlap:
    """Get the dict of overlaps between exons.

    Such dict has as keys a tuple of two connected nodes and as value the overlap between them:

    - Positive overlap means that they overlap that number of bases,
    - Zero that they occur next to each other
    - Negative that there is a gap in the transcriptome of that number of bases
        (one or multiple exons of length < kmer)

    :param splice_graph: returns: Note: the splice graph must have already the nodes written with
        coordinates, and the edges alredy entered too.

    Hypothesis: node2coords.values should only hold one value

    """
    logging.debug("\tComputing edge overlaps")

    # Init
    node2coords = nx.get_node_attributes(
        G=splice_graph,
        name='coordinates'
    )

    edge_overlaps = Edge2Overlap({edge: None for edge in splice_graph.edges()})

    for (node1, node2) in edge_overlaps.keys():
        node1_end = node2coords[node1][0][2]
        node2_start = node2coords[node2][0][1]

        # Overlap in bases, 0 means one next to the other, negative numbers a gap
        overlap = node1_end - node2_start
        edge_overlaps[(node1, node2)] = overlap

    return edge_overlaps


def build_splice_graph(bed6df: pd.DataFrame) -> SpliceGraph:
    """Build the splice_graph from a dataframe of bed6 records

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

    :param pd.DataFrame bed6df: Exon coordinates in BED6 format
    """
    logging.debug("Running build_splice_graph")
    # Initialize graph
    splice_graph = SpliceGraph()

    # Process nodes
    logging.debug("\tAdding nodes")
    splice_graph.add_nodes_from(bed6df["name"].tolist())
    nx.set_node_attributes(  # Add coordinates
        G=splice_graph,
        name="coordinates",
        values=bed6df_to_node2coordinates(bed6df)
    )

    # Process edges
    logging.debug("\tAdding edges")
    transcript2path = bed6df_to_path2node(bed6df)
    for path in transcript2path.values():
        splice_graph.add_path(path)
        nx.set_edge_attributes(
            G=splice_graph,
            name='overlaps',
            values=compute_edge_overlaps(splice_graph)
        )

    return splice_graph


def build_splice_graph_dict(
        bed3records: Iterable[Coordinate], args: dict) -> SpliceGraphDict:
    """Build the SpliceGraphDict from a bunch of BED3 records

    :param iterable bed3records: Iterable of tuples containing bed3 records
    :param dict args: args to be passed to the pipeline
    """
    logging.info("Building splice graph")

    # Process bed records
    bed6df_dict = bed3_records_to_bed6df_dict(bed3records)

    # Initialize pool of workers
    pool = mp.Pool(args["threads"])

    # Build graphs in parallel and merge results
    splice_graphs = pool.map(
        func=build_splice_graph,
        iterable=bed6df_dict.values(),
        chunksize=1000
    )
    pool.close()
    pool.join()
    splice_graph_dict = SpliceGraphDict(zip(bed6df_dict.keys(), splice_graphs))

    return splice_graph_dict
