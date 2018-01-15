#!/usr/bin/env python3

"""
exfi.io.gfa1_to_splice_graph.py: submodule to convert a gfa1 file into a splice graph
"""

import networkx as nx

from exfi.io.read_gfa1 import read_gfa1

def gfa1_to_splice_graph(handle):
    """(str) -> nx.DiGraph

    Read a GFA1 file and store the splice graph
    """

    # Read
    gfa1 = read_gfa1(handle)

    # Process
    coordinate_dict = gfa1["containments"]
    overlap_dict = gfa1["links"]

    splice_graph = nx.DiGraph()
    splice_graph.add_nodes_from(coordinate_dict.keys())
    splice_graph.add_edges_from(overlap_dict.keys())
    nx.set_node_attributes(
        G=splice_graph, name="coordinates", values=coordinate_dict
    )
    nx.set_edge_attributes(
        G=splice_graph, name="overlaps", values=overlap_dict
    )
