#!/usr/bin/env python3

"""tests.test_polish_overlaps.py: Tests for the exfi.polish submodule"""


from unittest import \
    TestCase, \
    main

import networkx as nx

from exfi.polish import \
    trim_start, \
    trim_end, \
    polish_overlaps, \
    polish_overlaps_dict

from tests.custom_assertions import \
    CustomAssertions

from tests.test_data import \
    SPLICE_GRAPH_EMPTY, SPLICE_GRAPH_SIMPLE, SPLICE_GRAPH_COMPLEX, \
    SPLICE_GRAPH_EMPTY_DICT, SPLICE_GRAPH_SIMPLE_DICT, SPLICE_GRAPH_COMPLEX_DICT, \
    TRANSCRIPTOME_EMPTY_DICT, TRANSCRIPTOME_SIMPLE_DICT, TRANSCRIPTOME_COMPLEX_DICT, \
    NODE2COORDS_COMPLEX_PART1, NODE2COORDS_COMPLEX_PART2, \
    OVERLAPS_COMPLEX, \
    SPLICE_GRAPH_1, SPLICE_GRAPH_2



# Test data
POLISHED_EMPTY = nx.DiGraph()
POLISHED_SIMPLE = nx.DiGraph()
POLISHED_SIMPLE.add_node("ENSDART00000161035.1:0-326")
nx.set_node_attributes(
    G=POLISHED_SIMPLE,
    name="coordinates",
    values={'ENSDART00000161035.1:0-326': (('ENSDART00000161035.1', 0, 326),)}
)


POLISHED_COMPLEX = nx.DiGraph()
POLISHED_COMPLEX.add_nodes_from(
    tuple(NODE2COORDS_COMPLEX_PART1.keys()) + tuple(NODE2COORDS_COMPLEX_PART2.keys())
)
nx.set_node_attributes(
    G=POLISHED_COMPLEX,
    name="coordinates",
    values={**NODE2COORDS_COMPLEX_PART1, **NODE2COORDS_COMPLEX_PART2}
)
POLISHED_COMPLEX.add_edges_from(OVERLAPS_COMPLEX.keys())
nx.set_edge_attributes(
    G=POLISHED_COMPLEX,
    name="overlaps",
    values=OVERLAPS_COMPLEX
)

POLISHED_EMPTY_DICT = dict()
POLISHED_SIMPLE_DICT = {"ENSDART00000161035.1": POLISHED_SIMPLE}
POLISHED_COMPLEX_DICT = {
    'ENSDART00000161035.1': SPLICE_GRAPH_1,
    'ENSDART00000165342.1': SPLICE_GRAPH_2
}

ARGS = {
    "threads": 1
}



class TestCoordAddLeft(TestCase):
    """Tests for exfi.polish.trim_start"""

    def test_empty(self):
        """exfi.polish.trim_start: empty case"""
        with self.assertRaises(IndexError):
            trim_start([], 1)

    def test_short(self):
        """exfi.polish.trim_start: short case"""
        with self.assertRaises(IndexError):
            trim_start([1], 1)

    def test_correct(self):
        """exfi.polish.trim_start: correct case"""
        actual = trim_start(("tr", 1, 2), 2)
        expected = ("tr", 3, 2)
        self.assertEqual(actual, expected)

    def test_too_long(self):
        """exfi.polish.trim_start: too long case"""
        actual = trim_start(("tr", 1, 2, 5, 5), 2)
        expected = ("tr", 3, 2, 5, 5)
        self.assertEqual(actual, expected)



class TestCoordAddRight(TestCase):
    """Tests for exfi.polish.trim_end"""

    def test_empty(self):
        """exfi.polish.trim_end: empty case"""
        with self.assertRaises(IndexError):
            trim_end([], 1)

    def test_short(self):
        """exfi.polish.trim_end: short case"""
        with self.assertRaises(IndexError):
            trim_end([1], 1)

    def test_correct(self):
        """exfi.polish.trim_end: correct case"""
        actual = trim_end(("tr", 1, 2), 2)
        expected = ("tr", 1, 0)
        self.assertEqual(actual, expected)

    def test_too_long(self):
        """exfi.polish.trim_end: too long case"""
        actual = trim_end(("tr", 1, 2, 5, 5), 2)
        expected = ("tr", 1, 0, 5, 5)
        self.assertEqual(actual, expected)



class TestPolishOverlaps(TestCase, CustomAssertions):
    """Tests for exfi.polish.polish_overlaps"""

    def test_empty(self):
        """exfi.polish.polish_overlaps: empty case"""
        actual = polish_overlaps(SPLICE_GRAPH_EMPTY, TRANSCRIPTOME_EMPTY_DICT)
        expected = POLISHED_EMPTY
        self.assertEqualSpliceGraphs(actual, expected)

    def test_simple(self):
        """exfi.polish.polish_overlaps: simple case"""
        actual = polish_overlaps(SPLICE_GRAPH_SIMPLE, TRANSCRIPTOME_SIMPLE_DICT)
        expected = POLISHED_SIMPLE
        self.assertEqualSpliceGraphs(actual, expected)

    def test_complex(self):
        """exfi.polish.polish_overlaps: complex case"""
        actual = polish_overlaps(SPLICE_GRAPH_COMPLEX, TRANSCRIPTOME_COMPLEX_DICT)
        expected = POLISHED_COMPLEX
        print(actual)
        print(expected)
        self.assertEqualSpliceGraphs(actual, expected)



class TestPolishOverlapsDict(TestCase, CustomAssertions):
    """Tests for exfi.polish.polish_overlaps_dict"""

    def test_empty(self):
        """exfi.polish.polish_overlaps_dict: empty case"""
        actual = polish_overlaps_dict(SPLICE_GRAPH_EMPTY_DICT, TRANSCRIPTOME_EMPTY_DICT, ARGS)
        expected = POLISHED_EMPTY_DICT
        self.assertEqualDictOfSpliceGraphs(actual, expected)

    def test_simple(self):
        """exfi.polish.polish_overlaps_dict: simple case"""
        actual = polish_overlaps_dict(SPLICE_GRAPH_SIMPLE_DICT, TRANSCRIPTOME_SIMPLE_DICT, ARGS)
        expected = POLISHED_SIMPLE_DICT
        self.assertEqualDictOfSpliceGraphs(actual, expected)

    def test_complex(self):
        """exfi.polish.polish_overlaps_dict: complex case"""
        actual = polish_overlaps_dict(SPLICE_GRAPH_COMPLEX_DICT, TRANSCRIPTOME_COMPLEX_DICT, ARGS)
        expected = POLISHED_COMPLEX_DICT
        self.assertEqualDictOfSpliceGraphs(actual, expected)



if __name__ == '__main__':
    main()
