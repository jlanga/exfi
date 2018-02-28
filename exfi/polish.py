#!/usr/bin/env python3

"""exfi.polish_overlaps: submodule to polish the overlaps between two exons"""

from typing import Iterable

import networkx as nx

from exfi.io import \
    _coordinate_to_str

from exfi.classes import Coordinate, SpliceGraph, SpliceGraphDict, FastaDict


def trim_end(coordinate: Coordinate, bases: int) -> Coordinate:
    """Trim bases to the end of the coordinate

    :param tuple coordinate: BED3 record
    :param int bases:  Number of bases to trim form end
    """
    return Coordinate(coordinate[0], coordinate[1], coordinate[2] - bases)



def trim_start(coordinate: Coordinate, bases: int) -> Coordinate:
    """Trim bases to the start of the coordinate

    :param tuple coordinate: BED3 record.
    :param int bases: Number of bases to trim from start.
    """
    return Coordinate(coordinate[0], coordinate[1] + bases, coordinate[2])



def trim_multiple_ends(iterable_coordinate: Iterable[Coordinate], bases: int) \
    -> Iterable[Coordinate]:
    """Trim bases at the end of all elements in iterable_coordinate

    :param tuple iterable_coordinate: iterable of bed3 records.
    :param int bases: number of bases to trim from end.
    """
    return tuple(trim_end(coordinate, bases) for coordinate in iterable_coordinate)



def trim_multiple_starts(iterable_coordinate: Iterable[Coordinate], bases: int) \
    -> Iterable[Coordinate]:
    """Trim bases at the start of all elements in iterable_coordinate

    :param tuple iterable_coordinate: iterable of bed3 records.
    :param int bases: Number of bases to trim from start.

    """
    return tuple(trim_start(coordinate, bases) for coordinate in iterable_coordinate)



def polish_splice_graph(splice_graph: SpliceGraph, fasta_dict: FastaDict) -> SpliceGraph:
    """Trim overlaps according to the AG/GT signal (AC/CT in the reverse strand)

    :param nx.DiGraph splice_graph: SpliceGraph to polish.
    :param dict fasta_dict: FastaDict of transcriptome.
    """

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
                new_node_u = _coordinate_to_str(trim_end(node_u_coord, overlap - index - 2))
                new_node_v = _coordinate_to_str(trim_start(node_v_coord, index + 2))

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



def polish_splice_graph_dict(
        splice_graph_dict: SpliceGraphDict, fasta_dict: FastaDict, args: dict) -> SpliceGraphDict:
    """Polish all overlaps in a splice graph dict

    :param dict splice_graph_dict: SpliceGraphDict to polish.
    :param dict fasta_dict: FastaDict of transcriptome.
    :param dict args: Dict of arguments for processing.

    args must at least be {"threads": 1}
    """

    import pathos.multiprocessing as mp

    # Initialize pool of workers
    pool = mp.Pool(args["threads"])

    # def polish_wrapper(splice_graph: SpliceGraph) -> SpliceGraph:
    #     """Export fasta_dict to function.
    #
    #     :param nx.DiGraph splice_graph: splice_graph to polish.
    #     """
    #     return polish_splice_graph(splice_graph, fasta_dict)

    # # Run
    # results = pool.starmap(
    #     func=polish_wrapper,
    #     iterable=splice_graph_dict.values(),
    #     chunksize=1000
    # )
    # pool.close()
    # pool.join()

    # Run
    results = pool.starmap(
        func=polish_splice_graph,
        iterable=(
            (splice_graph_value, {splice_graph_key: fasta_dict[splice_graph_key]})
            for splice_graph_key, splice_graph_value in splice_graph_dict.items()
        ),
        chunksize=1000
    )
    pool.close()
    pool.join()

    # Add results to splice_graph_dict
    for i, transcript in enumerate(splice_graph_dict.keys()):
        splice_graph_dict[transcript] = results[i]

    return splice_graph_dict
