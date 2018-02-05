#!/usr/bin/env python3

"""
Tests for exfi.build_splice_graph
"""


import unittest

import networkx as nx

from exfi.build_splice_graph import \
    _bed3_to_str, \
    bed3_records_to_bed6df_dict, \
    bed6df_to_path2node, \
    bed6df_to_node2coordinates, \
    compute_edge_overlaps, \
    build_splice_graph

from tests.custom_assertions import \
    CustomAssertions

from tests.test_data import \
    BED3RECORDS_EMPTY, BED3RECORDS_SIMPLE, BED3RECORDS_COMPLEX, \
    BED6DF_EMPTY, BED6DF_SIMPLE, BED6DF_COMPLEX, \
    BED6DF_DICT_EMPTY, BED6DF_DICT_SIMPLE, BED6DF_DICT_COMPLEX, \
    PATH_EMPTY, PATH_SIMPLE, PATH_COMPLEX, \
    NODE2COORDS_EMPTY, NODE2COORDS_SIMPLE, NODE2COORDS_COMPLEX, \
    OVERLAPS_EMPTY_DICT, OVERLAPS_SIMPLE_DICT, OVERLAPS_COMPLEX_DICT, \
    SPLICE_GRAPH_EMPTY_DICT, SPLICE_GRAPH_SIMPLE_DICT, SPLICE_GRAPH_COMPLEX_DICT

BED3_COLS = ['chrom', 'start', 'end']
BED6_COLS = ['chrom', 'start', 'end', 'name', 'score', 'strand']



def _prepare_overlaps(bed3_records):
    """Compute splicegraph prior the computation of overlaps"""
    sg_dict = dict()

    bed6df_dict = bed3_records_to_bed6df_dict(bed3_records)
    for transcript, bed6_df in bed6df_dict.items():
        splice_graph = nx.DiGraph()
        splice_graph.add_nodes_from(bed6_df["name"].tolist())
        node2coords = bed6df_to_node2coordinates(bed6_df)
        nx.set_node_attributes(G=splice_graph, name="coordinates", values=node2coords)
        transcript2path = bed6df_to_path2node(bed6_df)
        for path in transcript2path.values():
            splice_graph.add_path(path)
        sg_dict[transcript] = splice_graph
    return sg_dict



class TestBed3ToStr(unittest.TestCase):
    """Tests for _bed3_to_str"""

    def test_empty(self):
        """exfi.build_splice_graph._bed3_to_str: empty record"""
        with self.assertRaises(IndexError):
            _bed3_to_str([])

    def test_malformed1(self):
        """exfi.build_splice_graph._bed3_to_str: record of 2 elements"""
        with self.assertRaises(IndexError):
            _bed3_to_str((0, 1))

    def test_malformed2(self):
        """exfi.build_splice_graph._bed3_to_str: record of 4 elements"""
        with self.assertRaises(IndexError):
            _bed3_to_str((0, 1, 2, 3))

    def test_record(self):
        """exfi.build_splice_graph._bed3_to_str: correct record"""
        self.assertEqual(
            _bed3_to_str(("tr", 10, 15)),
            "tr:10-15"
        )



class TestBed3RecordsToBed6DFDict(unittest.TestCase, CustomAssertions):
    """Tests for bed3_records_to_bed6df_dict"""

    def test_empty_index(self):
        """exfi.build_splice_graph.bed3_records_to_bed6df_dict: empty exome"""
        actual = bed3_records_to_bed6df_dict(BED3RECORDS_EMPTY)
        expected = BED6DF_DICT_EMPTY
        self.assertEqualDictOfDF(actual, expected)

    def test_one_entry(self):
        """exfi.build_splice_graph.bed3_records_to_bed6df_dict: single exon"""
        actual = bed3_records_to_bed6df_dict(BED3RECORDS_SIMPLE)
        expected = BED6DF_DICT_SIMPLE
        self.assertEqualDictOfDF(actual, expected)

    def test_multiple(self):
        """exfi.build_splice_graph.bed3_records_to_bed6df_dict: multiple transcripts - multiple
        exons"""
        actual = bed3_records_to_bed6df_dict(BED3RECORDS_COMPLEX)
        expected = BED6DF_DICT_COMPLEX
        self.assertEqualDictOfDF(actual, expected)



class TestBed6DFToPath2Node(unittest.TestCase):
    """Tests for bed6df_to_path2node"""

    def test_empty(self):
        """exfi.build_splice_graph.bed6df_to_path2node: convert an empty exome to path"""
        actual = bed6df_to_path2node(BED6DF_EMPTY)
        expected = PATH_EMPTY
        self.assertEqual(actual, expected)

    def test_single(self):
        """exfi.build_splice_graph.bed6df_to_path2node: convert an single exon transcript to path"""
        actual = bed6df_to_path2node(BED6DF_SIMPLE)
        expected = PATH_SIMPLE
        self.assertEqual(actual, expected)

    def test_multiple(self):
        """exfi.build_splice_graph.bed6df_to_path2node: convert an single exon transcript to path"""
        actual = bed6df_to_path2node(BED6DF_COMPLEX)
        expected = PATH_COMPLEX
        self.assertEqual(actual, expected)



class TestBed6ToNode2Coord(unittest.TestCase):
    """Tests for bed6df_to_path2node"""

    def test_empty(self):
        """exfi.build_splice_graph.bed6df_to_node2coordinates: empty records"""
        actual = bed6df_to_node2coordinates(BED6DF_EMPTY)
        expected = NODE2COORDS_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """exfi.build_splice_graph.bed6df_to_node2coordinates: single node"""
        actual = bed6df_to_node2coordinates(BED6DF_SIMPLE)
        expected = NODE2COORDS_SIMPLE
        self.assertEqual(actual, expected)

    def test_complex(self):
        """exfi.build_splice_graph.bed6df_to_node2coordinates: complex case"""
        actual = bed6df_to_node2coordinates(BED6DF_COMPLEX)
        expected = NODE2COORDS_COMPLEX
        self.assertEqual(actual, expected)



class TestComputeEdgeOverlaps(unittest.TestCase):
    """Tests for compute_edge_overlaps"""

    def test_empty_exome(self):
        """exfi.build_splice_graph.compute_overlaps: compute the overlaps of an empty exome"""
        splice_graph_dict = _prepare_overlaps(BED3RECORDS_EMPTY)
        overlaps_dict = {
            transcript: compute_edge_overlaps(splice_graph)
            for transcript, splice_graph in splice_graph_dict.items()
        }
        self.assertEqual(overlaps_dict, OVERLAPS_EMPTY_DICT)

    def test_single_exon(self):
        """exfi.build_splice_graph.compute_overlaps: compute the overlaps of a single exon exome"""
        splice_graph_dict = _prepare_overlaps(BED3RECORDS_SIMPLE)
        overlaps_dict = {
            transcript: compute_edge_overlaps(splice_graph)
            for transcript, splice_graph in splice_graph_dict.items()
        }
        self.assertEqual(overlaps_dict, OVERLAPS_SIMPLE_DICT)

    def test_multiple_exons(self):
        """exfi.build_splice_graph.compute_overlaps: compute the overlaps of a simple exome"""
        splice_graph_dict = _prepare_overlaps(BED3RECORDS_COMPLEX)
        overlaps_dict = {
            transcript: compute_edge_overlaps(splice_graph)
            for transcript, splice_graph in splice_graph_dict.items()
        }
        self.assertEqual(overlaps_dict, OVERLAPS_COMPLEX_DICT)



class TestBuildSpliceGraph(unittest.TestCase, CustomAssertions):
    """Tests for build_splice_graph"""

    def test_empty(self):
        """exfi.build_splice_graph.build_splice_graph: compute the splice graph of an empty set
        of exons
        """
        actual = build_splice_graph(BED3RECORDS_EMPTY)
        expected = SPLICE_GRAPH_EMPTY_DICT
        self.assertEqualDictOfSpliceGraphs(actual, expected)

    def test_simple(self):
        """exfi.build_splice_graph.build_splice_graph: compute the splice graph of a singe exon"""
        actual = build_splice_graph(BED3RECORDS_SIMPLE)
        expected = SPLICE_GRAPH_SIMPLE_DICT
        self.assertEqualDictOfSpliceGraphs(actual, expected)

    def test_multiple(self):
        """exfi.build_splice_graph.build_splice_graph: compute the splice graph of a set of exons"""
        actual = build_splice_graph(BED3RECORDS_COMPLEX)
        expected = SPLICE_GRAPH_COMPLEX_DICT
        self.assertEqualDictOfSpliceGraphs(actual, expected)



if __name__ == '__main__':
    unittest.main()
