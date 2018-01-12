#!/usr/bin/env python3

"""
Tests for exfi.build_splice_graph
"""


import unittest

import networkx as nx

from exfi.build_splice_graph import \
    _bed3_to_str, \
    bed3_records_to_bed6df, \
    bed6df_to_path2node, \
    bed6df_to_node2coordinates, \
    compute_edge_overlaps, \
    build_splice_graph


from tests.test_data import \
    BED3RECORDS_EMPTY, BED3RECORDS_SIMPLE, BED3RECORDS_COMPLEX, \
    BED6DF_EMPTY, BED6DF_SIMPLE, BED6DF_COMPLEX, \
    PATH_EMPTY, PATH_SIMPLE, PATH_COMPLEX, \
    NODE2COORDS_EMPTY, NODE2COORDS_SIMPLE, NODE2COORDS_COMPLEX, \
    OVERLAPS_EMPTY, OVERLAPS_SIMPLE, OVERLAPS_COMPLEX, \
    SPLICE_GRAPH_EMPTY, SPLICE_GRAPH_SIMPLE, SPLICE_GRAPH_COMPLEX

BED3_COLS = ['chrom', 'start', 'end']
BED6_COLS = ['chrom', 'start', 'end', 'name', 'score', 'strand']


def _prepare_overlaps(bed3_records):
    """Compute splicegraph prior the computation of overlaps"""
    splice_graph = nx.DiGraph()
    bed6df = bed3_records_to_bed6df(bed3_records)
    splice_graph.add_nodes_from(bed6df["name"].tolist())
    node2coords = bed6df_to_node2coordinates(bed6df)
    nx.set_node_attributes(
        G=splice_graph,
        name="coordinates",
        values=node2coords
    )
    transcript2path = bed6df_to_path2node(bed6df)
    for path in transcript2path.values():
        splice_graph.add_path(path)
    return splice_graph


class TestBed3ToStr(unittest.TestCase):
    """Tests for _bed3_to_str"""
    def test_empty(self):
        """_bed3_to_str: empty record"""
        with self.assertRaises(IndexError):
            _bed3_to_str([])

    def test_malformed1(self):
        """_bed3_to_str: record of 2 elements"""
        with self.assertRaises(IndexError):
            _bed3_to_str((0, 1))

    def test_malformed2(self):
        """_bed3_to_str: record of 4 elements"""
        with self.assertRaises(IndexError):
            _bed3_to_str((0, 1, 2, 3))

    def test_record(self):
        """_bed3_to_str: correct record"""
        self.assertEqual(
            _bed3_to_str(("tr", 10, 15)),
            "tr:10-15"
        )


class TestBed3RecordsToBed6DF(unittest.TestCase):
    """Tests for bed3_records_to_bed6df"""

    def test_empty_index(self):
        """bed3_records_to_bed6df: empty exome"""
        actual = bed3_records_to_bed6df(BED3RECORDS_EMPTY)
        expected = BED6DF_EMPTY
        self.assertTrue(
            actual.equals(expected)
        )

    def test_one_entry(self):
        """bed3_records_to_bed6df: single exon"""
        actual = bed3_records_to_bed6df(BED3RECORDS_SIMPLE)
        expected = BED6DF_SIMPLE
        self.assertTrue(
            actual.equals(expected)
        )

    def test_multiple(self):
        """bed3_records_to_bed6df: multiple transcripts - multiple exons"""
        actual = bed3_records_to_bed6df(BED3RECORDS_COMPLEX)
        expected = BED6DF_COMPLEX
        self.assertTrue(
            actual.equals(expected)
        )


class TestBed6DFToPath2Node(unittest.TestCase):
    """Tests for bed6df_to_path2node"""

    def test_empty(self):
        """bed6df_to_path2node: convert an empty exome to path"""
        self.assertEqual(
            bed6df_to_path2node(BED6DF_EMPTY),
            PATH_EMPTY
        )


    def test_single(self):
        """bed6df_to_path2node: convert an single exon transcript to path"""
        self.assertEqual(
            bed6df_to_path2node(BED6DF_SIMPLE),
            PATH_SIMPLE
        )

    def test_multiple(self):
        """bed6df_to_path2node: convert an single exon transcript to path"""
        self.assertEqual(
            bed6df_to_path2node(BED6DF_COMPLEX),
            PATH_COMPLEX
        )


class TestBed6ToNode2Coord(unittest.TestCase):
    """Tests for bed6df_to_path2node"""
    def test_empty(self):
        """bed6df_to_node2coordinates: empty records"""
        self.assertEqual(
            NODE2COORDS_EMPTY,
            bed6df_to_node2coordinates(BED6DF_EMPTY)
        )

    def test_simple(self):
        """bed6df_to_node2coordinates: single node"""
        self.assertEqual(
            bed6df_to_node2coordinates(BED6DF_SIMPLE),
            NODE2COORDS_SIMPLE
        )

    def test_complex(self):
        """bed6df_to_node2coordinates: complex case"""
        self.assertEqual(
            bed6df_to_node2coordinates(BED6DF_COMPLEX),
            NODE2COORDS_COMPLEX
        )


class TestComputeEdgeOverlaps(unittest.TestCase):
    """Tests for compute_edge_overlaps"""
    def test_empty_exome(self):
        """compute_overlaps: compute the overlaps of an empty exome"""
        splice_graph = _prepare_overlaps(BED3RECORDS_EMPTY)
        overlaps = compute_edge_overlaps(splice_graph)
        self.assertEqual(overlaps, OVERLAPS_EMPTY)

    def test_single_exon(self):
        """compute_overlaps: compute the overlaps of a single exon exome"""
        splice_graph = _prepare_overlaps(BED3RECORDS_SIMPLE)
        overlaps = compute_edge_overlaps(splice_graph)
        self.assertEqual(overlaps, OVERLAPS_SIMPLE)

    def test_multiple_exons(self):
        """compute_overlaps: compute the overlaps of a simple exome"""
        splice_graph = _prepare_overlaps(BED3RECORDS_COMPLEX)
        overlaps = compute_edge_overlaps(splice_graph)
        self.assertEqual(overlaps, OVERLAPS_COMPLEX)


class TestBuildSpliceGraph(unittest.TestCase):
    """Tests for build_splice_graph"""

    def test_empty(self):
        """build_splice_graph: compute the splice graph of an empty set of exons"""
        actual = build_splice_graph(BED3RECORDS_EMPTY)
        self.assertTrue(
            nx.is_isomorphic(
                actual,
                SPLICE_GRAPH_EMPTY
            )
        )
        self.assertEqual(
            nx.get_node_attributes(G=actual, name="coordinates"),
            nx.get_node_attributes(G=SPLICE_GRAPH_EMPTY, name="coordinates"),
        )
        self.assertEqual(
            nx.get_edge_attributes(G=actual, name="overlaps"),
            nx.get_edge_attributes(G=SPLICE_GRAPH_EMPTY, name="overlaps"),
        )

    def test_simple(self):
        """build_splice_graph: compute the splice graph of a singe exon"""
        actual = build_splice_graph(BED3RECORDS_SIMPLE)
        self.assertTrue(nx.is_isomorphic(
            actual,
            SPLICE_GRAPH_SIMPLE
        ))
        self.assertEqual(
            nx.get_node_attributes(G=actual, name="coordinates"),
            nx.get_node_attributes(G=SPLICE_GRAPH_SIMPLE, name="coordinates"),
        )
        self.assertEqual(
            nx.get_edge_attributes(G=actual, name="overlaps"),
            nx.get_edge_attributes(G=SPLICE_GRAPH_SIMPLE, name="overlaps"),
        )

    def test_multiple(self):
        """build_splice_graph: compute the splice graph of a set of exons"""
        actual = build_splice_graph(BED3RECORDS_COMPLEX)
        self.assertTrue(nx.is_isomorphic(
            actual,
            SPLICE_GRAPH_COMPLEX
        ))
        self.assertEqual(
            nx.get_node_attributes(G=actual, name="coordinates"),
            nx.get_node_attributes(G=SPLICE_GRAPH_COMPLEX, name="coordinates"),
        )
        self.assertEqual(
            nx.get_edge_attributes(G=actual, name="overlaps"),
            nx.get_edge_attributes(G=SPLICE_GRAPH_COMPLEX, name="overlaps"),
        )


if __name__ == '__main__':
    unittest.main()
