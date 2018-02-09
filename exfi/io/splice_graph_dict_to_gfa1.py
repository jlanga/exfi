#!/usr/bin/env python3

"""
exfi.io.splice_graph_to_gfa1.py: functions to convert a splice graph into a GFA1 file
"""

import logging

from itertools import chain

import networkx as nx
import pandas as pd

from natsort import \
    natsorted

from exfi.build_splice_graph_dict import \
    bed6df_to_path2node



def _get_node2coord(splice_graph):
    """Get node coordinates"""
    return nx.get_node_attributes(G=splice_graph, name="coordinates")



def _get_edge2overlap(splice_graph):
    """Get edge overlap"""
    return nx.get_edge_attributes(G=splice_graph, name="overlaps")



def _compute_segments(splice_graph_dict, transcriptome_dict):
    """(nx.DiGraph, dict) -> list

    Compute the segment lines:
    S node_id sequence length
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



def _compute_links(splice_graph_dict: dict):
    """(dict of nx.DiGraph) -> list

    Compute the link lines:
    L start orientation end orientation overlap
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



def _compute_containments(splice_graph_dict):
    """(dict of nx.DiGraph) -> list

    Compute the containment lines (w.r.t. transcriptome):
    C container orientation contained orientation position overlap
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



def _compute_paths(splice_graph_dict):
    """(dict of nx.DiGraph) -> lists

    Compute the paths in the splice graph:
    P transcript_id [node1, ..., nodeN]
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



def splice_graph_dict_to_gfa1(splice_graph_dict, transcriptome_dict, filename):
    """(dict_of_seqrecords, dict_of_seqrecords, str) -> None

    Write splice graph to filename in GFA 1 format
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
