#!/usr/bin/env python3

"""
Tests for exfi.io.gfa1_to_splice_graph
"""

from unittest import TestCase, main

import networkx as nx

from exfi.io.gfa1_to_splice_graph import gfa1_to_splice_graph

from tests.test_data import \
    SPLICE_GRAPH_EMPTY, SPLICE_GRAPH_SIMPLE, SPLICE_GRAPH_COMPLEX, \
    GFA_EMPTY_FN, GFA_SIMPLE_FN, GFA_COMPLEX_FN



class TestGFA1ToSpliceGrah(TestCase):
    """Tests for exfi.io.gfa1_to_splice_graph.gfa1_to_splice_graph"""

    def test_empty(self):
        """exfi.io.gfa1_to_splice_graph.gfa1_to_splice_graph: empty case"""
        actual = gfa1_to_splice_graph(GFA_EMPTY_FN)
        self.assertTrue(nx.is_isomorphic(
            actual,
            SPLICE_GRAPH_EMPTY
        ))
        self.assertEqual(
            nx.get_node_attributes(G=actual, name="coordinates"),
            nx.get_node_attributes(G=SPLICE_GRAPH_EMPTY, name="coordinates"),
        )
        self.assertEqual(
            nx.get_edge_attributes(G=actual, name="overlaps"),
            nx.get_edge_attributes(G=SPLICE_GRAPH_EMPTY, name="overlaps"),
        )


    def test_simple(self):
        """exfi.io.gfa1_to_splice_graph.gfa1_to_splice_graph: simple case"""
        actual = gfa1_to_splice_graph(GFA_SIMPLE_FN)
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


    def test_complex(self):
        """exfi.io.gfa1_to_splice_graph.gfa1_to_splice_graph: complex case"""
        actual = gfa1_to_splice_graph(GFA_COMPLEX_FN)
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
    main()
