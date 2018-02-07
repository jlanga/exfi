#!/usr/bin/env python3

"""
Tests for exfi.collapse_splice_graph
"""

from unittest import TestCase, main

from exfi.collapse_splice_graph_dict import \
    _compute_seq2node, \
    _compute_old2new, \
    _compute_new_node2coord, \
    _compute_new_link2overlap, \
    collapse_splice_graph_dict

from tests.custom_assertions import \
    CustomAssertions

from tests.test_data import \
    NODE2COORDS_EMPTY, NODE2COORDS_SIMPLE, NODE2COORDS_COMPLEX, \
    TRANSCRIPTOME_EMPTY_DICT, TRANSCRIPTOME_SIMPLE_DICT, TRANSCRIPTOME_COMPLEX_DICT, \
    SPLICE_GRAPH_EMPTY_DICT, SPLICE_GRAPH_SIMPLE_DICT, SPLICE_GRAPH_COMPLEX_DICT, \
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
        actual = _compute_seq2node(NODE2COORDS_EMPTY, TRANSCRIPTOME_EMPTY_DICT)
        expected = SEQ2NODE_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """exfi.collapse_splice_graph._compute_seq2node: simple case"""
        actual = _compute_seq2node(NODE2COORDS_SIMPLE, TRANSCRIPTOME_SIMPLE_DICT)
        expected = SEQ2NODE_SIMPLE
        self.assertEqual(actual, expected)

    def test_complex(self):
        """exfi.collapse_splice_graph._compute_seq2node: complex case"""
        actual = _compute_seq2node(NODE2COORDS_COMPLEX, TRANSCRIPTOME_COMPLEX_DICT)
        expected = SEQ2NODE_COMPLEX
        self.assertEqual(actual, expected)



class TestComputeOld2New(TestCase):
    """Tests for exfi.collapse_splice_graph._compute_old2new"""

    def test_empty(self):
        """exfi.collapse_splice_graph._compute_old2new: empty case"""
        actual = _compute_old2new(SEQ2NODE_EMPTY)
        expected = OLD2NEW_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """exfi.collapse_splice_graph._compute_seq2node: simple case"""
        actual = _compute_old2new(SEQ2NODE_SIMPLE)
        expected = OLD2NEW_SIMPLE
        self.assertEqual(actual, expected)

    def test_complex(self):
        """exfi.collapse_splice_graph._compute_seq2node: complex case"""
        actual = _compute_old2new(SEQ2NODE_COMPLEX)
        expected = OLD2NEW_COMPLEX
        self.assertEqual(actual, expected)



class TestComputeNewNode2Coord(TestCase):
    """Tests for exfi.collapse_splice_graph._compute_new_node2coord"""

    def test_empty(self):
        """exfi.collapse_splice_graph._compute_new_node2coord: empty case"""
        actual = _compute_new_node2coord(OLD2NEW_EMPTY, NODE2COORDS_EMPTY)
        expected = NEW_NODE2COORD_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """exfi.collapse_splice_graph._compute_new_node2coord: simple case"""
        actual = _compute_new_node2coord(OLD2NEW_SIMPLE, NODE2COORDS_SIMPLE)
        expected = NEW_NODE2COORD_SIMPLE
        self.assertEqual(actual, expected)

    def test_complex(self):
        """exfi.collapse_splice_graph._compute_new_node2coord: complex case"""
        actual = _compute_new_node2coord(OLD2NEW_COMPLEX, NODE2COORDS_COMPLEX)
        expected = NEW_NODE2COORD_COMPLEX
        self.assertEqual(actual, expected)



class TestComputeNewLink2Overlap(TestCase):
    """Tests for exfi.collapse_splice_graph._compute_new_link2overlap"""

    def test_empty(self):
        """exfi.collapse_splice_graph._compute_new_link2overlap: empty case"""
        actual = _compute_new_link2overlap(OLD2NEW_EMPTY, LINK2OVERLAP_EMPTY)
        expected = NEW_LINK2OVERLAP_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """exfi.collapse_splice_graph._compute_new_link2overlap: simple case"""
        actual = _compute_new_link2overlap(OLD2NEW_SIMPLE, LINK2OVERLAP_SIMPLE)
        expected = NEW_LINK2OVERLAP_SIMPLE
        self.assertEqual(actual, expected)

    def test_complex(self):
        """exfi.collapse_splice_graph._compute_new_link2overlap: complex case"""
        actual = _compute_new_link2overlap(OLD2NEW_COMPLEX, LINK2OVERLAP_COMPLEX)
        expected = NEW_LINK2OVERLAP_COMPLEX
        self.assertEqual(actual, expected)




class TestCollapseSpliceGraph(TestCase, CustomAssertions):
    """Tests for exfi.collapse_splice_graph_dict.collapse_splice_graph_dict"""

    def test_empty(self):
        """exfi.collapse_splice_graph.collapse_splice_graph: empty case"""
        actual = collapse_splice_graph_dict(SPLICE_GRAPH_EMPTY_DICT, TRANSCRIPTOME_EMPTY_DICT)
        expected = COLLAPSED_EMPTY
        self.assertEqualSpliceGraphs(actual, expected)

    def test_simple(self):
        """exfi.collapse_splice_graph.collapse_splice_graph: simple case"""
        actual = collapse_splice_graph_dict(SPLICE_GRAPH_SIMPLE_DICT, TRANSCRIPTOME_SIMPLE_DICT)
        expected = COLLAPSED_SIMPLE
        self.assertEqualSpliceGraphs(actual, expected)

    def test_complex(self):
        """exfi.collapse_splice_graph.collapse_splice_graph: complex case"""
        actual = collapse_splice_graph_dict(SPLICE_GRAPH_COMPLEX_DICT, TRANSCRIPTOME_COMPLEX_DICT)
        expected = COLLAPSED_COMPLEX
        self.assertEqualSpliceGraphs(actual, expected)



if __name__ == "__main__":
    main()
