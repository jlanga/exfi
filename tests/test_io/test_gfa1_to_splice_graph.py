#!/usr/bin/env python3

"""
Tests for exfi.io.gfa1_to_splice_graph
"""

from unittest import TestCase, main

from exfi.io.gfa1_to_splice_graph_dict import gfa1_to_splice_graph_dict

from tests.custom_assertions import \
    CustomAssertions

from tests.test_data import \
    SPLICE_GRAPH_EMPTY_DICT, SPLICE_GRAPH_SIMPLE_DICT, SPLICE_GRAPH_COMPLEX_DICT, \
    GFA_EMPTY_FN, GFA_SIMPLE_FN, GFA_COMPLEX_FN



class TestGFA1ToSpliceGrah(TestCase, CustomAssertions):
    """Tests for exfi.io.gfa1_to_splice_graph.gfa1_to_splice_graph"""

    def test_empty(self):
        """exfi.io.gfa1_to_splice_graph.gfa1_to_splice_graph: empty case"""
        actual = gfa1_to_splice_graph_dict(GFA_EMPTY_FN)
        expected = SPLICE_GRAPH_EMPTY_DICT
        self.assertEqualDictOfSpliceGraphs(actual, expected)


    def test_simple(self):
        """exfi.io.gfa1_to_splice_graph.gfa1_to_splice_graph: simple case"""
        actual = gfa1_to_splice_graph_dict(GFA_SIMPLE_FN)
        expected = SPLICE_GRAPH_SIMPLE_DICT
        self.assertEqualDictOfSpliceGraphs(actual, expected)


    def test_complex(self):
        """exfi.io.gfa1_to_splice_graph.gfa1_to_splice_graph: complex case"""
        actual = gfa1_to_splice_graph_dict(GFA_COMPLEX_FN)
        expected = SPLICE_GRAPH_COMPLEX_DICT
        self.assertEqualDictOfSpliceGraphs(actual, expected)


if __name__ == '__main__':
    main()
