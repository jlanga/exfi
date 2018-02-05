#!/usr/bin/env python3

"""
Tests for exfi.io.gfa1_to_splice_graph
"""

from unittest import TestCase, main

from exfi.io.gfa1_to_splice_graph import gfa1_to_splice_graph

from tests.custom_assertions import \
    CustomAssertions

from tests.test_data import \
    SPLICE_GRAPH_EMPTY, SPLICE_GRAPH_SIMPLE, SPLICE_GRAPH_COMPLEX, \
    GFA_EMPTY_FN, GFA_SIMPLE_FN, GFA_COMPLEX_FN



class TestGFA1ToSpliceGrah(TestCase, CustomAssertions):
    """Tests for exfi.io.gfa1_to_splice_graph.gfa1_to_splice_graph"""

    def test_empty(self):
        """exfi.io.gfa1_to_splice_graph.gfa1_to_splice_graph: empty case"""
        actual = gfa1_to_splice_graph(GFA_EMPTY_FN)
        self.assertEqualSpliceGraphs(
            actual,
            SPLICE_GRAPH_EMPTY
        )


    def test_simple(self):
        """exfi.io.gfa1_to_splice_graph.gfa1_to_splice_graph: simple case"""
        actual = gfa1_to_splice_graph(GFA_SIMPLE_FN)
        self.assertEqualSpliceGraphs(
            actual,
            SPLICE_GRAPH_SIMPLE
        )


    def test_complex(self):
        """exfi.io.gfa1_to_splice_graph.gfa1_to_splice_graph: complex case"""
        actual = gfa1_to_splice_graph(GFA_COMPLEX_FN)
        self.assertEqualSpliceGraphs(
            actual,
            SPLICE_GRAPH_COMPLEX
        )


if __name__ == '__main__':
    main()
