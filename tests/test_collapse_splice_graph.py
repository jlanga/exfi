#!/usr/bin/env python3

"""
Tests for exfi.collapse_splice_graph
"""

from unittest import TestCase, main

from exfi.collapse_splice_graph import \
    _compute_seq2node, \
    _compute_old2new, \
    _compute_new_node2coord, \
    _compute_new_link2overlap, \
    collapse_splice_graph

from tests.custom_assertions import \
    CustomAssertions

from tests.test_data import \
    NODE2COORDS_EMPTY, NODE2COORDS_SIMPLE, NODE2COORDS_COMPLEX, \
    TRANSCRIPTOME_EMPTY_DICT, TRANSCRIPTOME_SIMPLE_DICT, TRANSCRIPTOME_COMPLEX_DICT, \
    SPLICE_GRAPH_EMPTY, SPLICE_GRAPH_SIMPLE, SPLICE_GRAPH_COMPLEX, \
    SEQ2NODE_EMPTY, SEQ2NODE_SIMPLE, SEQ2NODE_COMPLEX, \
    OLD2NEW_EMPTY, OLD2NEW_SIMPLE, OLD2NEW_COMPLEX, \
    NEW_NODE2COORD_EMPTY, NEW_NODE2COORD_SIMPLE, NEW_NODE2COORD_COMPLEX, \
    LINK2OVERLAP_EMPTY, LINK2OVERLAP_SIMPLE, LINK2OVERLAP_COMPLEX, \
    NEW_LINK2OVERLAP_EMPTY, NEW_LINK2OVERLAP_SIMPLE, NEW_LINK2OVERLAP_COMPLEX, \
    COLLAPSED_EMPTY, COLLAPSED_SIMPLE, COLLAPSED_COMPLEX



class TestComputeSeq2Node(TestCase):
    """Tests for exfi.collapse_splice_graph._compute_seq2node"""

    def test_empty(self):
        """exfi.collapse_splice_graph._compute_seq2node: empty case"""
        self.assertEqual(
            _compute_seq2node(NODE2COORDS_EMPTY, TRANSCRIPTOME_EMPTY_DICT),
            SEQ2NODE_EMPTY
        )

    def test_simple(self):
        """exfi.collapse_splice_graph._compute_seq2node: simple case"""
        self.assertEqual(
            _compute_seq2node(NODE2COORDS_SIMPLE, TRANSCRIPTOME_SIMPLE_DICT),
            SEQ2NODE_SIMPLE
        )

    def test_complex(self):
        """exfi.collapse_splice_graph._compute_seq2node: complex case"""
        self.assertEqual(
            _compute_seq2node(NODE2COORDS_COMPLEX, TRANSCRIPTOME_COMPLEX_DICT),
            SEQ2NODE_COMPLEX
        )



class TestComputeOld2New(TestCase):
    """Tests for exfi.collapse_splice_graph._compute_old2new"""

    def test_empty(self):
        """exfi.collapse_splice_graph._compute_old2new: empty case"""
        self.assertEqual(
            _compute_old2new(SEQ2NODE_EMPTY),
            OLD2NEW_EMPTY
        )

    def test_simple(self):
        """exfi.collapse_splice_graph._compute_seq2node: simple case"""
        self.assertEqual(
            _compute_old2new(SEQ2NODE_SIMPLE),
            OLD2NEW_SIMPLE
        )

    def test_complex(self):
        """exfi.collapse_splice_graph._compute_seq2node: complex case"""
        self.assertEqual(
            _compute_old2new(SEQ2NODE_COMPLEX),
            OLD2NEW_COMPLEX
        )



class TestComputeNewNode2Coord(TestCase):
    """Tests for exfi.collapse_splice_graph._compute_new_node2coord"""

    def test_empty(self):
        """exfi.collapse_splice_graph._compute_new_node2coord: empty case"""
        self.assertEqual(
            _compute_new_node2coord(OLD2NEW_EMPTY, NODE2COORDS_EMPTY),
            NEW_NODE2COORD_EMPTY
        )

    def test_simple(self):
        """exfi.collapse_splice_graph._compute_new_node2coord: simple case"""
        self.assertEqual(
            _compute_new_node2coord(OLD2NEW_SIMPLE, NODE2COORDS_SIMPLE),
            NEW_NODE2COORD_SIMPLE
        )

    def test_complex(self):
        """exfi.collapse_splice_graph._compute_new_node2coord: complex case"""
        self.assertEqual(
            _compute_new_node2coord(OLD2NEW_COMPLEX, NODE2COORDS_COMPLEX),
            NEW_NODE2COORD_COMPLEX
        )



class TestComputeNewLink2Overlap(TestCase):
    """Tests for exfi.collapse_splice_graph._compute_new_link2overlap"""

    def test_empty(self):
        """exfi.collapse_splice_graph._compute_new_link2overlap: empty case"""
        self.assertEqual(
            _compute_new_link2overlap(OLD2NEW_EMPTY, LINK2OVERLAP_EMPTY),
            NEW_LINK2OVERLAP_EMPTY
        )

    def test_simple(self):
        """exfi.collapse_splice_graph._compute_new_link2overlap: simple case"""
        self.assertEqual(
            _compute_new_link2overlap(OLD2NEW_SIMPLE, LINK2OVERLAP_SIMPLE),
            NEW_LINK2OVERLAP_SIMPLE
        )

    def test_complex(self):
        """exfi.collapse_splice_graph._compute_new_link2overlap: complex case"""
        self.assertEqual(
            _compute_new_link2overlap(OLD2NEW_COMPLEX, LINK2OVERLAP_COMPLEX),
            NEW_LINK2OVERLAP_COMPLEX
        )




class TestCollapseSpliceGraph(TestCase, CustomAssertions):
    """Tests for exfi.collapse_splice_graph.collapse_splice_graph"""

    def test_empty(self):
        """exfi.collapse_splice_graph.collapse_splice_graph: empty case"""
        self.assertEqualSpliceGraphs(
            collapse_splice_graph(SPLICE_GRAPH_EMPTY, TRANSCRIPTOME_EMPTY_DICT),
            COLLAPSED_EMPTY
        )

    def test_simple(self):
        """exfi.collapse_splice_graph.collapse_splice_graph: simple case"""
        self.assertEqualSpliceGraphs(
            collapse_splice_graph(SPLICE_GRAPH_SIMPLE, TRANSCRIPTOME_SIMPLE_DICT),
            COLLAPSED_SIMPLE
        )

    def test_complex(self):
        """exfi.collapse_splice_graph.collapse_splice_graph: complex case"""
        self.assertEqualSpliceGraphs(
            collapse_splice_graph(SPLICE_GRAPH_COMPLEX, TRANSCRIPTOME_COMPLEX_DICT),
            COLLAPSED_COMPLEX
        )



if __name__ == "__main__":
    main()
