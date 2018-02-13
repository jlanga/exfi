#!/usr/bin/env python3

"""exfi.polish_overlaps: submodule to polish the overlaps between two exons"""

import networkx as nx

from exfi.build_splice_graph_dict import \
    _bed3_to_str


def trim_end(coordinate: tuple, bases: int) -> tuple:
    """Trim bases to the end of the coordinate"""
    coordinate = list(coordinate)
    coordinate[2] -= bases
    return tuple(coordinate)



def trim_start(coordinate: tuple, bases: int) -> tuple:
    """Trim bases to the start of the coordinate"""
    coordinate = list(coordinate)
    coordinate[1] += bases
    return tuple(coordinate)



def trim_multiple_ends(iterable_coordinate: tuple, bases: int) -> tuple:
    """Trim bases at the end of all elements in iterable_coordinate"""
    return tuple(trim_end(coordinate, bases) for coordinate in iterable_coordinate)



def trim_multiple_starts(iterable_coordinate: tuple, bases: int) -> tuple:
    """Trim bases at the start of all elements in iterable_coordinate"""
    return tuple(trim_start(coordinate, bases) for coordinate in iterable_coordinate)



def polish_overlaps(splice_graph, fasta_dict):
    """Trim overlaps according to the AG/GT signal (AC/CT in the reverse strand)"""

    node2coordinates = nx.get_node_attributes(G=splice_graph, name="coordinates")
    edge2overlap = nx.get_edge_attributes(G=splice_graph, name="overlaps")
    node_mapping = {node: node for node in splice_graph.nodes()}  # old: new

    for (node_u, node_v), overlap in edge2overlap.items():

        if overlap >= 4:

            # Get one of the coordinates (there should be one)
            node_u_coord = node2coordinates[node_u][0]
            node_v_coord = node2coordinates[node_v][0]

            # Get overlapping thing
            overlap_seq = fasta_dict[node_u_coord[0]][node_v_coord[1]:node_u_coord[2]]

            # When there is an overlap,
            # Exon structure should be EXON...AG - GT...intron...AG - GT...exon
            if "AGGT" in overlap_seq:

                index = overlap_seq.rfind("AGGT")

                # rename both transcripts
                ## u: delete overlap untul AG
                ## v: delete overlap until GT
                new_node_u = _bed3_to_str(trim_end(node_u_coord, overlap - index - 2))
                new_node_v = _bed3_to_str(trim_start(node_v_coord, index + 2))

                # Update old -> new renaming
                node_mapping[node_u] = new_node_u
                node_mapping[node_v] = new_node_v

                # change coordinate dict values
                ## u
                node2coordinates[node_u] = trim_multiple_ends(
                    iterable_coordinate=node2coordinates[node_u], bases=overlap - index - 2
                )
                ## v
                node2coordinates[node_v] = trim_multiple_starts(
                    iterable_coordinate=node2coordinates[node_v], bases=index + 2
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
