#!/usr/bin/env python3

"""
Submodule to split a huge splice graph into its constituent components
"""

import logging

import networkx as nx
from natsort import natsorted


def split_into_components(splice_graph: nx.DiGraph) -> dict:
    '''Convert a single splice graph into a dict of splice graphs, where
    - keys are transcript_ids
    - values are the splice graph of that transcript
    '''
    logging.info("\tSplitting directed graph into directed components")
    # Compute connected components
    logging.info("\t\tComputing undirected components")
    undirected_components = nx.connected_component_subgraphs(G=splice_graph.to_undirected())

    component_dict = {}

    logging.info("\t\tComputing directed components")
    for undirected_component in undirected_components:

        # Get the transcript_id of the component
        nodes = tuple(x for x in undirected_component.nodes())
        a_node = nodes[0]
        transcript = undirected_component.node[a_node]["coordinates"][0][0]
        logging.info("\t\t\tProcessing component {transcript}".format(transcript=transcript))
        # Get node data as is
        node2coord = nx.get_node_attributes(
            G=undirected_component,
            name="coordinates"
        )

        # Get edge data. Be careful because each edge is twice: one in each direction
        # Use natsorted to get the correct direction
        edge2overlap = {}
        for node_u, node_v in undirected_component.edges():
            node_u, node_v = natsorted([node_u, node_v])
            edge2overlap[(node_u, node_v)] = splice_graph[node_u][node_v]["overlaps"]

        # Re-create directed graph
        # Nodes
        directed_component = nx.DiGraph()
        directed_component.add_nodes_from(node2coord.keys())
        nx.set_node_attributes(
            G=directed_component,
            name="coordinates",
            values=node2coord
        )
        # Edges
        directed_component.add_edges_from(edge2overlap.keys())
        nx.set_edge_attributes(
            G=directed_component,
            name="overlaps",
            values=edge2overlap
        )

        # Store directed component in its position
        component_dict[transcript] = directed_component

    return component_dict
