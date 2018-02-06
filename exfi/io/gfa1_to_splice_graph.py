#!/usr/bin/env python3

"""
exfi.io.gfa1_to_splice_graph.py: submodule to convert a gfa1 file into a splice graph
"""

import logging

import networkx as nx

from exfi.io.read_gfa1 import \
    read_gfa1


def _split_node2coord(node2coord: dict, node2transcript: dict) -> dict:
    """Split the big node2coord dict into its subcompoenents (transcripts)"""
    splitted_node2coord = {key: dict() for key in set(node2transcript.values())}
    for node, coordinates in node2coord.items():
        transcript = node2transcript[node]
        if node not in splitted_node2coord[transcript]:
            splitted_node2coord[transcript][node] = tuple()
        splitted_node2coord[transcript][node] += coordinates
    return splitted_node2coord



def _split_edge2overlap(edge2overlap: dict, node2transcript: dict) -> dict:
    """Split the big edge2overlap dict into subcomponents (transcripts)"""
    splitted_edge2overlap = {transcript: dict() for transcript in set(node2transcript.values())}
    for edge, overlap in edge2overlap.items():
        transcript = node2transcript[edge[0]]
        splitted_edge2overlap[transcript][edge] = overlap
    return splitted_edge2overlap



def gfa1_to_splice_graph(handle):
    """(str) -> nx.DiGraph

    Read a GFA1 file and store the splice graph
    """
    logging.info("Converting gfa1 %s to splice graph", handle)

    # Read and process
    logging.info("\tReading and processing GFA file %s", handle)
    gfa1 = read_gfa1(handle)
    node2coord = gfa1["containments"]
    edge2overlap = gfa1["links"]
    transcript2nodes = gfa1["paths"]

    # Revert path2nodes
    node2transcript = {
        value: key
        for key, values in transcript2nodes.items()
        for value in values
    }

    # Split node2coord
    transcript2node2coord = _split_node2coord(node2coord, node2transcript)
    transcript2edge2overlap = _split_edge2overlap(edge2overlap, node2transcript)

    # Initialize
    splice_graph_dict = {transcript: None for transcript in transcript2nodes}

    # process
    for transcript in splice_graph_dict.keys():
        splice_graph = nx.DiGraph()

        node2coord = transcript2node2coord[transcript]
        splice_graph. add_nodes_from(node2coord.keys())
        nx.set_node_attributes(
            G=splice_graph, name="coordinates", values=node2coord
        )

        edge2overlap = transcript2edge2overlap[transcript]
        splice_graph.add_edges_from(edge2overlap.keys())
        nx.set_edge_attributes(
            G=splice_graph, name="overlaps", values=edge2overlap
        )

        splice_graph_dict[transcript] = splice_graph


    return splice_graph_dict
