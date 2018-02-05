#!/usr/bin/env python3

"""exfi.io.join_components: Submodule to convert a dict of splice graphs into a
single splice graph"""

import logging

import networkx as nx



def join_components(dict_of_components: dict) -> nx.DiGraph:
    """Merge all splice graphs in dict_of_components into a single splice_graph"""

    logging.info("\tJoining multiple splice graphs into one")

    # Join everything into a splice_graph
    joint = nx.DiGraph()

    # Nodes
    logging.info("\t\tProcessing nodes")
    node2coordinate = {
        node: coordinate
        for subgraph in dict_of_components.values()
        for node, coordinate in nx.get_node_attributes(G=subgraph, name="coordinates").items()
    }

    joint.add_nodes_from(node2coordinate.keys())
    nx.set_node_attributes(G=joint, name="coordinates", values=node2coordinate)

    # Edges
    logging.info("\t\tProcessing edges")
    edge2overlap = {
        edge: overlap
        for subgraph in dict_of_components.values()
        for edge, overlap in nx.get_edge_attributes(G=subgraph, name="overlaps").items()
    }
    joint.add_edges_from(edge2overlap.keys())
    nx.set_edge_attributes(G=joint, name="overlaps", values=edge2overlap)

    return joint
