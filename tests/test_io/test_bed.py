#!/usr/bin/ebv python3

"""tests.io.bed.py: tests for exfi.io.bed.py"""

import unittest

from exfi.io.bed import \
    bed3_to_bed4, \
    bed4_to_node2coordinates, \
    bed4_to_path2nodes, \
    bed4_to_node2sequence, \
    bed4_to_edge2overlap

from tests.io.bed import \
    BED3_EMPTY, BED3_SIMPLE, BED3_COMPLEX, \
    BED4_EMPTY, BED4_SIMPLE, BED4_COMPLEX, \
    NODE2COORDINATES_EMPTY, NODE2COORDINATES_SIMPLE, NODE2COORDINATES_COMPLEX, \
    PATH2NODES_EMPTY, PATH2NODES_SIMPLE, PATH2NODES_COMPLEX, \
    NODE2SEQUENCE_EMPTY, NODE2SEQUENCE_SIMPLE, NODE2SEQUENCE_COMPLEX, \
    EDGE2OVERLAP_EMPTY, EDGE2OVERLAP_SIMPLE, EDGE2OVERLAP_COMPLEX

from tests.io.transcriptome_dicts import \
    TRANSCRIPTOME_EMPTY_DICT, TRANSCRIPTOME_SIMPLE_DICT, \
    TRANSCRIPTOME_COMPLEX_DICT

class TestBed3ToBed4(unittest.TestCase):
    """Tests for exfi.io.bed.bed3_to_bed4"""

    def test_empty(self):
        """exfi.io.bed.bed3_to_bed4: empty case"""
        observed = bed3_to_bed4(BED3_EMPTY)
        self.assertTrue(BED4_EMPTY.equals(observed))

    def test_simple(self):
        """exfi.io.bed.bed3_to_bed4: simple case"""
        observed = bed3_to_bed4(BED3_SIMPLE)
        self.assertTrue(observed.equals(BED4_SIMPLE))

    def test_complex(self):
        """exfi.io.bed.bed3_to_bed4: complex case"""
        observed = bed3_to_bed4(BED3_COMPLEX)
        self.assertTrue(observed.equals(BED4_COMPLEX))


class TestBed4ToNode2Coordinates(unittest.TestCase):
    """Tests for exfi.io.bed.bed4_to_node2coordinates"""
    def test_empty(self):
        """exfi.io.bed.bed4_to_node2coordinates: empty case"""
        observed = bed4_to_node2coordinates(BED4_EMPTY)
        self.assertTrue(observed.equals(NODE2COORDINATES_EMPTY))

    def test_simple(self):
        """exfi.io.bed.bed4_to_node2coordinates: simple case"""
        observed = bed4_to_node2coordinates(BED4_SIMPLE)
        self.assertTrue(observed.equals(NODE2COORDINATES_SIMPLE))

    def test_complex(self):
        """exfi.io.bed.bed4_to_node2coordinates: complex case"""
        observed = bed4_to_node2coordinates(BED4_COMPLEX)
        self.assertTrue(observed.equals(NODE2COORDINATES_COMPLEX))


class TestBed4ToPath2Nodes(unittest.TestCase):
    """Tests for exfi.io.bed.bed4_to_path2nodes"""
    def test_empty(self):
        """exfi.io.bed.bed4_to_path2nodes: empty case"""
        observed = bed4_to_path2nodes(BED4_EMPTY)
        self.assertEqual(observed, PATH2NODES_EMPTY)

    def test_simple(self):
        """exfi.io.bed.bed4_to_path2nodes: simple case"""
        observed = bed4_to_path2nodes(BED4_SIMPLE)
        self.assertEqual(observed, PATH2NODES_SIMPLE)

    def test_complex(self):
        """exfi.io.bed.bed4_to_path2nodes: complex case"""
        observed = bed4_to_path2nodes(BED4_COMPLEX)
        print("Observed:\n", observed)
        print("Expected:\n", PATH2NODES_COMPLEX)
        self.assertEqual(observed, PATH2NODES_COMPLEX)


class TestBed4ToNode2Sequence(unittest.TestCase):
    """Tests for exfi.io.bed.bed4_to_node2sequence"""

    def test_empty(self):
        """exfi.io.bed.bed4_to_node2sequence: empty case"""
        observed = bed4_to_node2sequence(BED4_EMPTY, TRANSCRIPTOME_EMPTY_DICT)
        self.assertTrue(observed.equals(NODE2SEQUENCE_EMPTY))

    def test_simple(self):
        """exfi.io.bed.bed4_to_node2sequence: simple case"""
        observed = bed4_to_node2sequence(BED4_SIMPLE, TRANSCRIPTOME_SIMPLE_DICT)
        self.assertTrue(observed.equals(NODE2SEQUENCE_SIMPLE))

    def test_complex(self):
        """exfi.io.bed.bed4_to_node2sequence: complex case"""
        observed = bed4_to_node2sequence(
            BED4_COMPLEX, TRANSCRIPTOME_COMPLEX_DICT
        )
        self.assertTrue(observed.equals(NODE2SEQUENCE_COMPLEX))


class TestBed4ToEdge2Overlap(unittest.TestCase):
    """Tests for exfi.io.bed.bed4_to_edge2overlap"""

    def test_empty(self):
        """exfi.io.bed.bed4_to_edge2overlap: empty case"""
        observed = bed4_to_edge2overlap(BED4_EMPTY)
        self.assertTrue(observed.equals(EDGE2OVERLAP_EMPTY))

    def test_simple(self):
        """exfi.io.bed.bed4_to_edge2overlap: simple case"""
        observed = bed4_to_edge2overlap(BED4_SIMPLE)
        print("Observed:\n", observed, observed.dtypes)
        print("Expected:\n", EDGE2OVERLAP_SIMPLE, EDGE2OVERLAP_SIMPLE.dtypes)
        self.assertTrue(observed.equals(EDGE2OVERLAP_SIMPLE))

    def test_complex(self):
        """exfi.io.bed.bed4_to_edge2overlap: complex case"""
        observed = bed4_to_edge2overlap(BED4_COMPLEX)
        print("Observed:\n", observed, observed.dtypes)
        print("Expected:\n", EDGE2OVERLAP_COMPLEX, EDGE2OVERLAP_COMPLEX.dtypes)
        self.assertTrue(observed.equals(EDGE2OVERLAP_COMPLEX))


if __name__ == '__main__':
    unittest.main()
