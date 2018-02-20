#!/usr/bin/env python3.6

"""Classes for exfi:

- FastaDict
- Node2Coordinates
- Edge2Overlap
- SpliceGraph
- SpliceGraphDict
"""



from typing import \
    Tuple, Dict, Iterable

from Bio.SeqIO.FastaIO import \
    SimpleFastaParser

import networkx as nx

from exfi.correct import correct_splice_graph_dict
from exfi.collapse import collapse_splice_graph_dict
from exfi.polish import polish_splice_graph_dict
from exfi.io.splice_graph_dict_to_gfa1 import splice_graph_dict_to_gfa1
from exfi.io.gfa1_to_splice_graph_dict import gfa1_to_splice_graph_dict
from exfi.build_splice_graph_dict import build_splice_graph_dict



class FastaDict(Dict[str, str]):
    """Class for fasta dictionaries.

    FastaDict is basically a dictionary where keys are sequence indentifiers and values nucleotide
    sequences as str.
    """
    # pylint: disable=too-few-public-methods

    @staticmethod
    def build_from_fasta(filename: str) -> None:
        """Read Fasta file"""
        fasta_dict = FastaDict()
        with open(filename, "r") as handle:
            for key, value in SimpleFastaParser(handle):
                fasta_dict[key] = value
        return fasta_dict



class Coordinate(Tuple[str, int, int]):
    """Store coordinates as tuples of the shape str, int, int"""
    seqid = None
    start = None
    end = None

    def __init__(self, seqid: str, start: int, end: int):
        super(Coordinate, self).__init__()
        if isinstance(seqid, str) and isinstance(start, int) and isinstance(end, int):
            self.seqid = seqid
            self.start = start
            self.end = end
        else:
            type_seqid = type(seqid)
            type_start = type(start)
            type_end = type(end)
            raise TypeError(
                f"seqid must be str: {type_seqid}, "
                f"start must be int: {type_start}, and end must be"
                f" int: {type_end}."
            )

    def __str__(self) -> None:
        """Print to : - notation

        ("chr1", 100, 150) -> "chr1:100-150"
        """
        print(f"{self.seqid}:{self.start}-{self.end}")

    def exted_bases_at_start(self, bases):
        """Extend bases at start

        ("chr1", 100, 150) + 5 -> ("chr1", 95, 150)
        """
        self.start -= bases

    def extend_bases_at_end(self, bases):
        """Extend bases at end"

        ("chr1", 100, 150) + 5 -> ("chr", 100, 155)
        """
        self.end += bases



class Node2Coordinates(Dict[str, Tuple[Coordinate]]):
    """Node2Coordinates is a dict where keys are node ids and values is a tuple of tuples containing
    the exon coordinates with respect to the transcriptome in form (transcript, start, end)."""
    # pylint: disable=too-few-public-methods
    pass



class Edge2Overlap(Dict[Tuple[str, str], int]):
    """Edge2Overlaps is a dict where keys are tuples of two node identifiers and the value is an int
    """
    # pylint: disable=too-few-public-methods
    pass



class SpliceGraph(nx.DiGraph):
    """Class to work with single splice graphs.

    A SpliceGraph is basically a nx.DiGraph whose:
        - nodes are exon identifiers,
        - links are nodes that are contiguous,
        - nodes have an attribute called "coordinates", which is a Node2Coordinates object.
        - edges have an attribute called "overlap", which is a Edge2Overlap object.
    """

    def __init__(self) -> None:
        super(SpliceGraph).__init__()
        self = nx.DiGraph()

    def set_node2coordinates(self, node2coordinates: Node2Coordinates = None) -> None:
        """Use nx.set_node_attributes"""
        if isinstance(node2coordinates, Node2Coordinates):
            nx.set_node_attributes(G=self.splice_graph, name="coordinates", values=node2coordinates)
        else:
            raise TypeError("node2coordinates is not Node2Coordinates: ", type(node2coordinates))

    def get_node2coordinates(self) -> Node2Coordinates:
        """Use nx.get_node_attributes and convert to Node2Coordinate."""
        return Node2Coordinates(nx.get_node_attributes(G=self.splice_graph, name="cooridnates"))

    def add_edges(self, edges: Iterable[Tuple[str, str]]) -> None:
        """Add edges via nx.add_edges_from"""
        self.splice_graph.add_edges_from(edges)

    def set_edge2overlap(self, edge2overlap: Edge2Overlap = None) -> None:
        """Set edge to overlap values"""
        if isinstance(edge2overlap, Edge2Overlap):
            nx.set_edge_attributes(G=self.splice_graph, name="overlap", values=edge2overlap)
        else:
            raise TypeError("edge2overlap is not Edge2Overlap: ", type(edge2overlap))

    def get_edge2overlap(self) -> Edge2Overlap:
        """Get edge2overlap data from SpliceGraph."""
        return Edge2Overlap(nx.get_edge_attributes(G=self.splice_graph, name="overlaps"))



class SpliceGraphDict(Dict[str, SpliceGraph]):
    """Class to work with splice graphs in dict format "name": SpliceGraph"""

    def add_from_iterables(
            self, names: Iterable[str], splice_graphs: Iterable[SpliceGraph]) -> None:
        """Adds to SpliceGraphDict the splice_graphs_i with name names_i"""
        for name, splice_graph in zip(names, splice_graphs):
            self[name] = splice_graph

    def polish(self, fasta_dict: FastaDict, args: dict) -> None:
        """Polish overlaps according to the signal AGGT is found"""
        self = polish_splice_graph_dict(splice_graph_dict=self, fasta_dict=fasta_dict, args=args)

    def correct(self, args: dict) -> None:
        """Use abyss-sealer to correct/merge exons"""
        self = correct_splice_graph_dict(splice_graph_dict=self, args=args)

    def collapse(self, transcriptome_dict: FastaDict) -> None:
        """Create a new splice graph by merging all exons by sequence identity"""
        self = collapse_splice_graph_dict(
            splice_graph_dict=self, transcriptome_dict=transcriptome_dict)

    def write_to_gfa1(self, transcriptome_dict: FastaDict, filename: str) -> None:
        """Write SpliceGraphDict to GFA1 file"""
        splice_graph_dict_to_gfa1(
            splice_graph_dict=self, transcriptome_dict=transcriptome_dict, filename=filename)

    @staticmethod
    def load_from_bed3_records(
            bed3_records: Iterable[Coordinate], args: dict) -> None:
        """Build the SpliceGraphDict from BED3 records"""
        return SpliceGraphDict(build_splice_graph_dict(bed3records=bed3_records, args=args))

    @staticmethod
    def load_from_gfa1_file(filename) -> None:
        """Build SpliceGraphDict from a GFA1 file"""
        return SpliceGraphDict(gfa1_to_splice_graph_dict(handle=filename))



# class Path2Node(dict):

# class Bed3Records

# class Bed3Df

# class Bed6Df
