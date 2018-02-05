#!/usr/bin/env python3

"""Tests for exfi.io.split_into_components"""

from unittest import TestCase, main

from exfi.io.join_components import join_components
from exfi.io.split_into_components import split_into_components

from tests.custom_assertions import \
    CustomAssertions

from tests.test_data import \
    SPLICE_GRAPH_EMPTY, SPLICE_GRAPH_SIMPLE, SPLICE_GRAPH_COMPLEX



class TestSplitIntoComponents(TestCase, CustomAssertions):
    """Tests for exfi.io.split_into_components.split_into_components

     split_into_components(splice_graph: nx.DiGraph) -> dict
    """

    def test_empty(self):
        """exfi.io.split_into_components.split_into_components: empty case"""
        actual = join_components(split_into_components(SPLICE_GRAPH_EMPTY))
        self.assertEqualSpliceGraphs(actual, SPLICE_GRAPH_EMPTY)

    def test_simple(self):
        """exfi.io.split_into_components.split_into_components: simple case"""
        actual = join_components(split_into_components(SPLICE_GRAPH_SIMPLE))
        self.assertEqualSpliceGraphs(actual, SPLICE_GRAPH_SIMPLE)

    def test_complex(self):
        """exfi.io.split_into_components.split_into_components: complex case"""
        actual = join_components(split_into_components(SPLICE_GRAPH_COMPLEX))
        self.assertEqualSpliceGraphs(actual, SPLICE_GRAPH_COMPLEX)



if __name__ == '__main__':
    main()
