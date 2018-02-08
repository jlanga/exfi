#!/usr/bin/env python3

"""exfi.polish_overlaps: submodule to polish the overlaps between two exons"""

import networkx as nx

from exfi.build_splice_graph_dict import \
    _bed3_to_str


def coord_add_left(coord, bases):
    """Add bases to the start coordinate"""
    coord = list(coord)
    coord[1] += bases
    return tuple(coord)



def coord_add_right(coord, bases):
    """Add bases to the end coordinate"""
    coord = list(coord)
    coord[2] += bases
    return tuple(coord)



def polish_overlaps(splice_graph, fasta_dict):
    """Trim overlaps according to the AG/GT signal (AC/CT in the reverse strand)"""

    node2coordinates = nx.get_node_attributes(G=splice_graph, name="coordinates")
    edge2overlap = nx.get_edge_attributes(G=splice_graph, name="overlaps")
    node_mapping = {node: node for node in splice_graph.nodes()}  # old: new

    for (node_u, node_v), overlap in edge2overlap.items():

        if overlap >= 3:

            # Get coordinates
            node_u_coord = node2coordinates[node_u][0]
            node_v_coord = node2coordinates[node_v][0]

            # Get overlapping thing
            overlap_seq = fasta_dict[node_u_coord[0]][node_v_coord[1]:node_u_coord[2]]
            difference = overlap_seq.find("GT") - overlap_seq.find("AG")

            # if AG.*GT:
            if difference > 0:

                # rename both transcripts (don't insert and delete)
                ## u
                new_node_u = _bed3_to_str(coord_add_right(node_u_coord, difference))
                new_node_v = _bed3_to_str(coord_add_left(node_v_coord, difference))

                # Update old -> new
                node_mapping[node_u] = new_node_u
                node_mapping[node_v] = new_node_v

                # change coordinate dict values
                ## u
                node2coordinates[node_u] = tuple(
                    coord_add_right(coordinate, difference)
                    for coordinate in node2coordinates[node_u]
                )
                ## v
                node2coordinates[node_v] = tuple(
                    coord_add_left(coordinate, difference)
                    for coordinate in  node2coordinates[node_v]
                )

                # change overlap dict values
                edge2overlap[(node_u, node_v)] = 0

            # else:
                # merge nodes into one
                # Leave as it is?

    # rename nodes in graph
    splice_graph = nx.relabel_nodes(G=splice_graph, mapping=node_mapping)

    # assign attributes
    ## Nodes
    nx.set_node_attributes(
        G=splice_graph,
        name="coordinates",
        values={
            node_mapping[node]: coordinates
            for node, coordinates in node2coordinates.items()
        }
    )
    ## Edges
    nx.set_edge_attributes(
        G=splice_graph,
        name="overlaps",
        values={
            (node_mapping[node_u], node_mapping[node_v]): overlap
            for (node_u, node_v), overlap in edge2overlap.items()
        }
    )

    return splice_graph



def polish_overlaps_dict(splice_graph_dict: dict, fasta_dict: dict, args: dict) -> dict:
    """Polish all overlaps in a splice graph dict"""

    import pathos.multiprocessing as mp

    # Initialize pool of workers
    pool = mp.Pool(args["threads"])

    def polish_wrapper(splice_graph):
        """Export all fasta to function"""
        return polish_overlaps(splice_graph, fasta_dict)

    # Run
    results = pool.map(
        func=polish_wrapper,
        iterable=splice_graph_dict.values(),
        chunksize=1000
    )

    # Add results to splice_graph_dict
    for i, transcript in enumerate(splice_graph_dict.keys()):
        splice_graph_dict[transcript] = results[i]

    return splice_graph_dict
