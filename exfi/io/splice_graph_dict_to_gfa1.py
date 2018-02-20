#!/usr/bin/env python3

"""exfi.io.splice_graph_to_gfa1.py: functions to convert a splice graph into a GFA1 file"""

import logging
from typing import Generator, Dict, Tuple

from itertools import chain

import networkx as nx
import pandas as pd

from natsort import \
    natsorted

from exfi.build_splice_graph_dict import \
    bed6df_to_path2node


def _get_node2coord(splice_graph: nx.DiGraph) -> Dict[str, Tuple[Tuple[str, int, int]]]:
    """Get node coordinates

    :param nx.DiGraph splice_graph: DiGraph from which to extract node coordinates.
    """
    return nx.get_node_attributes(G=splice_graph, name="coordinates")


def _get_edge2overlap(splice_graph: nx.DiGraph) -> Dict[Tuple[str, str], int]:
    """Get edge2overlap

    :param nx.DiGraph splice_graph: DiGraph from which to extract edge overlaps.
    """
    return nx.get_edge_attributes(G=splice_graph, name="overlaps")


def _set_node2coord(
        splice_graph: nx.DiGraph, node2coord: Dict[str, Tuple[Tuple[str, int, int]]]) -> None:
    """Set node to coordinates data in splice_graph

    :param nx.DiGraph splice_graph: DiGraph where to store the data.
    :param dict node2coord: node coordinates to be stored.
    """
    nx.set_node_attributes(G=splice_graph, name="coordinates", values=node2coord)


def _set_edge2overlap(splice_graph: nx.DiGraph, edge2overlap: Dict[Tuple[str, str], int]) -> None:
    """Set edge to overlap data in splice_graph

    :param nx.DiGraph splice_graph: DiGraph where to store data
    :param dict edge2overlap: edge to overlap data to store
    """
    nx.set_edge_attributes(G=splice_graph, name="overlaps", values=edge2overlap)


def _compute_segments(
        splice_graph_dict: Dict[str, nx.DiGraph], transcriptome_dict: Dict[str, str]) -> \
        Generator[str, None, None]:
    """Compute the segment lines: S node_id sequence length

    :param dict splice_graph_dict: dict of the shape {component_id: nx.DiGraph}
    :param dict transcriptome_dict: dict with the transcriptome FASTA
    """
    logging.info("\tComputing segments")
    for _, splice_graph in natsorted(splice_graph_dict.items()):
        node2coords = _get_node2coord(splice_graph)
        for node_id, coordinates in natsorted(node2coords.items()):
            logging.debug("\t\tProcessing node %s", node_id)
            coordinate = coordinates[0]
            transcript_id, start, end = coordinate
            sequence = str(transcriptome_dict[transcript_id][start:end])
            yield "S\t{node}\t{sequence}\tLN:i:{length}\n".format(
                node=node_id,
                sequence=sequence,
                length=len(sequence)
            )


def _compute_links(splice_graph_dict: Dict[str, nx.DiGraph]) -> Generator[str, None, None]:
    """Compute the link lines: L start orientation end orientation overlap

    :param dict splice_graph_dict: dict of DiGraphs with the splice graph
    """
    logging.info("\tComputing links")

    for _, splice_graph in natsorted(splice_graph_dict.items()):

        # Edges
        edge2overlap = _get_edge2overlap(splice_graph)

        for (node1, node2), overlap in natsorted(edge2overlap.items()):
            logging.debug("\t\tProcesssing edge (%s, %s)", node1, node2)
            # is an overlap or a gap
            if overlap >= 0:
                overlap = "{}M".format(overlap)
            else:
                overlap = "{}G".format(-overlap)

            yield "L\t{node1}\t{orientation1}\t{node2}\t{orientation2}\t{overlap}\n".format(
                node1=node1, orientation1="+",
                node2=node2, orientation2="+",
                overlap=overlap
            )


def _compute_containments(splice_graph_dict: Dict[str, nx.DiGraph]) -> Generator[str, None, None]:
    """Compute the containment lines (w.r.t. transcriptome)

    C container orientation contained orientation position overlap

    :param dict splice_graph_dict: dict of DiGraph representing the splice graph.
    """
    # Extract from the graph necessary data
    logging.info("\tComputing containments")

    for _, splice_graph in natsorted(splice_graph_dict.items()):
        node2coordinates = _get_node2coord(splice_graph)

        for node, coordinates in natsorted(node2coordinates.items()):
            for (transcript_id, start, end) in coordinates:
                logging.debug("\t\tProcessing %s - %s:%s-%s", node, transcript_id, start, end)
                cigar = str(int(end) - int(start)) + "M"
                yield "C\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                    transcript_id, "+", node, "+", start, cigar
                )


def _compute_paths(splice_graph_dict: Dict[str, nx.DiGraph]):
    """Compute the paths in the splice graph: P transcript_id [node1, ..., nodeN]

    :param splice_graph_dict: Dict of splice graphs.
    """
    logging.info("\tComputing paths")

    for _, splice_graph in natsorted(splice_graph_dict.items()):

        # Make bed6df
        node2coords = _get_node2coord(splice_graph)

        bed6_records = (
            (seqid, start, end, name, ".", "+")
            for name, coordinates in node2coords.items()
            for (seqid, start, end) in coordinates
        )

        bed6df = pd.DataFrame(
            data=bed6_records,
            columns=["chrom", "start", "end", "name", "score", "strand"]
        )

        path2nodes = bed6df_to_path2node(bed6df)
        for transcript_id, path in natsorted(path2nodes.items()):
            yield "P\t{transcript_id}\t{path}\n".format(
                transcript_id=transcript_id,
                path=",".join([node + "+" for node in path])
            )


def splice_graph_dict_to_gfa1(
        splice_graph_dict: Dict[str, nx.DiGraph],
        transcriptome_dict: Dict[str, str],
        filename: str) \
        -> None:
    """Write splice graph to filename in GFA 1 format

    :param splice_graph_dict: Dict of Splice Graphs.
    :param transcriptome_dict: Dict of Sequences.
    :param filename: Ouptut filename.
    """
    logging.info("Writing splice graph to GFA1 file %s", filename)
    header = ["H\tVN:Z:1.0\n"]
    segments = _compute_segments(splice_graph_dict, transcriptome_dict)
    links = _compute_links(splice_graph_dict)
    containments = _compute_containments(splice_graph_dict)
    paths = _compute_paths(splice_graph_dict)
    with open(filename, "w") as gfa1_out:
        gfa1_out.writelines(chain(
            header, segments, links, containments, paths
        ))
