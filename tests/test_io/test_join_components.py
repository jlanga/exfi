#!/usr/bin/env python3

"""Tests for exfi.io.join_components submodule"""

from unittest import TestCase, main

from exfi.io.join_components import join_components
from exfi.io.split_into_components import split_into_components
from exfi.io.gfa1_to_splice_graph import gfa1_to_splice_graph

from tests.auxiliary_functions import CustomAssertions

from tests.test_data import \
    SPLICE_GRAPH_EMPTY, SPLICE_GRAPH_SIMPLE, SPLICE_GRAPH_COMPLEX

class TestJoinComponents(TestCase, CustomAssertions):
    """Tests for exfi.io.join_components.join_components"""

    def test_empty(self):
        """exfi.io.join_components.join_components: empty case"""
        actual = join_components(split_into_components(SPLICE_GRAPH_EMPTY))
        self.assertEqualSpliceGraphs(actual, SPLICE_GRAPH_EMPTY)


    def test_simple(self):
        """exfi.io.join_components.join_components: simple case"""
        actual = join_components(split_into_components(SPLICE_GRAPH_SIMPLE))
        self.assertEqualSpliceGraphs(actual, SPLICE_GRAPH_SIMPLE)

    def test_complex(self):
        """exfi.io.join_components.join_components: complex case"""
        actual = join_components(split_into_components(SPLICE_GRAPH_EMPTY))
        self.assertEqualSpliceGraphs(actual, SPLICE_GRAPH_SIMPLE)





if __name__ == '__main__':
    main()
