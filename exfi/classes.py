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

class FastaDict(Dict[str, str]):
    """Class for fasta dictionaries.

    FastaDict is basically a dictionary where keys are sequence indentifiers and values nucleotide
    sequences as str.
    """
    # pylint: disable=too-few-public-methods

    def build_from_fasta(self, filename: str) -> None:
        """Read Fasta file"""
        with open(filename, "r") as handle:
            for key, value in SimpleFastaParser(handle):
                self[key] = value

    # def __str__(self) -> None:
    #     """Print the entire dictionary"""
    #     print("{")
    #     for key, value in self.items():
    #         print(f"{key}: {value},")
    #     print("}")



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
        self.start += bases

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



class SpliceGrapDict(Dict[str, SpliceGraph]):
    """Class to work with splice"""

    def add_from_iterables(
            self, names: Iterable[str], splice_graphs: Iterable[SpliceGraph]) -> None:
        """Adds to SpliceGraphDict the splice_graphs_i with name names_i"""
        for name, splice_graph in zip(names, splice_graphs):
            self[name] = splice_graph

    def polish(self, args: dict) -> None:
        """Polish overlaps according to the signal AGGT is found"""
        pass

    def correct(self, args: dict) -> None:
        """Use abyss-sealer to correct/merge exons"""
        pass

    def collapse(self, args: dict) -> None:
        """Create a new splice graph by merging all exons by sequence identity"""
        pass

    def write_to_gfa1(self, filename: str) -> None:
        """Write SpliceGraphDict to GFA1 file"""
        pass

    def load_from_bed3_records(self, bed3_records: Iterable[Tuple[str, int, int]]) -> None:
        """Build the SpliceGraphDict from BED3 records"""

    def load_from_gfa1_file(self, filename) -> None:
        """Build SpliceGraphDict from a GFA1 file"""

    def write_exons(self, filename) -> None:
        """Convert SpliceGraphDict to exons FASTA."""

    def write_gapped(self, filename) -> None:
        """Convert SpliceGraphDict to gapped transcripts in FASTA format."""



# class Path2Node(dict):

# class Bed3Records

# class Bed3Df

# class Bed6Df
