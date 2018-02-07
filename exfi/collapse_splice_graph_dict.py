#!/usr/bin/env python3

"""
exfi.collapse_splice_graph: merge nodes with the exact same sequence
"""

import logging

import networkx as nx

def _compute_seq2node(node2coord, transcriptome_dict):
    """(dict, dict) -> dict

    From
    - a dict node2coord = {node_id: tuples of coordinates}, and
    - a dict transcriptome_dict = {transcript_id: str},
    build the dict
    - node2seq = {node_id: str}
    """
    # Get the node -> sequence
    logging.debug("\t_compute_seq2node")
    seq2node = {}
    for node_id, coordinates in node2coord.items():
        seqid, start, end = coordinates[0]  # Just take the first
        sequence = transcriptome_dict[seqid][start:end]
        if sequence not in seq2node:
            seq2node[sequence] = ()
        seq2node[sequence] += (node_id, )
    return seq2node



def _compute_old2new(seq2node):
    """(dict) -> dict

    Compute the dict of old identifiers to new
    """
    logging.debug("\t_compute_old2new")
    old2new = {}
    for i, old_nodes in enumerate(seq2node.values()):
        new_node = "exon_{exon_number:08d}".format(exon_number=i)
        for old_node in old_nodes:
            old2new[old_node] = new_node
    return old2new



def _compute_new_node2coord(old2new, node2coord):
    """(dict, dict) -> dict

    Recompute the node to coordinate dict
    """
    logging.debug("\t_compute_new_node2coord")
    # Compute the new set coordinates of each node
    new_node2coord = {}
    for old_id, new_id in old2new.items():
        if new_id not in new_node2coord:
            new_node2coord[new_id] = ()
        new_node2coord[new_id] += node2coord[old_id]
    return new_node2coord


def _compute_new_link2overlap(old2new, link2overlap):
    """(dict, dict) -> dict

    Recompute the link2overlaps dict accordint to the new node_ids
    """
    logging.debug("\t_compute_new_link2overlap")
    # Compute the new set of edges and overlaps
    new_link2overlap = {}
    for (node_from, node_to), overlap in link2overlap.items():
        new_from = old2new[node_from]
        new_to = old2new[node_to]
        new_link2overlap[(new_from, new_to)] = overlap
    return new_link2overlap


def _merge_node2coords(splice_graph_dict):
    """Take all node2coord from every graph and merge into a single dict"""
    logging.debug("\t_merge_node2coords")
    node2coord_big = {}
    for splice_graph in splice_graph_dict.values():
        # Get node and edge data
        node2coord = nx.get_node_attributes(G=splice_graph, name="coordinates")

        for node, coordinate in node2coord.items():
            if node not in node2coord_big:
                node2coord_big[node] = tuple()
            node2coord_big[node] += coordinate

    return node2coord_big



def _merge_link2overlap(splice_graph_dict):
    """Take all link2overlap from every graph and merge into a single dict"""

    logging.debug("\t_merge_link2overlap")

    link2overlap_big = {}

    for splice_graph in splice_graph_dict.values():
        link2overlap = nx.get_edge_attributes(G=splice_graph, name="overlaps")
        for edge, overlap in link2overlap.items():
            link2overlap_big[edge] = overlap

    return link2overlap_big



def collapse_splice_graph_dict(splice_graph_dict, transcriptome_dict):
    """(dict nx.DiGraph, dict) -> nx.DiGraph

    Collapse nodes by sequence identity
    """
    logging.info("Collapsing graph by sequence")

    node2coord = _merge_node2coords(splice_graph_dict)
    link2overlap = _merge_link2overlap(splice_graph_dict)

    seq2node = _compute_seq2node(node2coord, transcriptome_dict)
    old2new = _compute_old2new(seq2node)
    del seq2node

    new_node2coord = _compute_new_node2coord(old2new, node2coord)
    new_link2overlap = _compute_new_link2overlap(old2new, link2overlap)

    # Build graph
    collapsed_graph = nx.DiGraph()
    collapsed_graph.add_nodes_from(new_node2coord.keys())
    collapsed_graph.add_edges_from(new_link2overlap.keys())
    nx.set_node_attributes(G=collapsed_graph, name="coordinates", values=new_node2coord)
    nx.set_edge_attributes(G=collapsed_graph, name="overlaps", values=new_link2overlap)

    return collapsed_graph
