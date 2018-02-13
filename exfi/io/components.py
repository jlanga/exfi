#!/usr/bin/env python3

"""exfi.io.join_components: Submodule to convert a dict of splice graphs into a
single splice graph"""

import logging

import networkx as nx

from natsort import natsorted

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

    node2coord_big = nx.get_node_attributes(G=splice_graph, name="coordinates")
    edge2overlap_big = nx.get_edge_attributes(G=splice_graph, name="overlaps")

    logging.info("\t\tComputing directed components")
    for undirected_component in undirected_components:

        # Get the transcript_id of the component
        logging.info("\t\t\tGetting component info")
        nodes = tuple(x for x in undirected_component.nodes())
        a_node = nodes[0]
        transcript = undirected_component.node[a_node]["coordinates"][0][0]
        logging.info("\t\t\tProcessing component %s", transcript)

        logging.info("\t\t\t\tGetting node2coord")
        node2coord = {
            node: node2coord_big[node]
            for node in undirected_component.nodes()
        }

        logging.info("\t\t\t\tGetting edge2overlap")
        edges = {
            tuple(natsorted([node_u, node_v]))
            for node_u, node_v in undirected_component.edges()
        }
        edge2overlap = {
            edge: edge2overlap_big[edge]
            for edge in edges
        }

        # Re-create directed graph
        # Nodes
        logging.info("\t\t\t\tCreating empty graph")
        directed_component = nx.DiGraph()

        logging.info("\t\t\t\tAdding nodes")
        directed_component.add_nodes_from(node2coord.keys())
        nx.set_node_attributes(
            G=directed_component,
            name="coordinates",
            values=node2coord
        )
        # Edges
        logging.info("\t\t\t\tAdding edges")
        directed_component.add_edges_from(edge2overlap.keys())
        nx.set_edge_attributes(
            G=directed_component,
            name="overlaps",
            values=edge2overlap
        )

        # Store directed component in its position
        component_dict[transcript] = directed_component

    return component_dict
