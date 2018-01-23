#!/usr/bin/env python3

"""
exfi.io.gfa1_to_splice_graph.py: submodule to convert a gfa1 file into a splice graph
"""

import logging

import networkx as nx

from exfi.io.read_gfa1 import read_gfa1

def gfa1_to_splice_graph(handle):
    """(str) -> nx.DiGraph

    Read a GFA1 file and store the splice graph
    """
    logging.info("Converting gfa1 {gfa} to splice graph".format(gfa=handle))

    # Read
    gfa1 = read_gfa1(handle)

    # Process
    logging.info("\tProcessing coordinates")
    coordinate_dict = gfa1["containments"]
    logging.info("\tProcessing links")
    overlap_dict = gfa1["links"]

    # Build
    logging.info("\tBuilding splice graph")
    splice_graph = nx.DiGraph()
    splice_graph.add_nodes_from(coordinate_dict.keys())
    splice_graph.add_edges_from(overlap_dict.keys())
    nx.set_node_attributes(
        G=splice_graph, name="coordinates", values=coordinate_dict
    )
    nx.set_edge_attributes(
        G=splice_graph, name="overlaps", values=overlap_dict
    )

    return splice_graph
