#!/usr/bin/env python3

"""tests.test_io.test_read_gfa.py: tests for exfi.io.read_gfa.py"""


from unittest import TestCase, main

from exfi.io.read_gfa import read_gfa1

from tests.io.gfa1 import \
    HEADER, \
    SEGMENTS_EMPTY, SEGMENTS_SIMPLE, SEGMENTS_COMPLEX, \
    SEGMENTS_COMPLEX_SOFT, SEGMENTS_COMPLEX_HARD, \
    LINKS_EMPTY, LINKS_SIMPLE, LINKS_COMPLEX, \
    CONTAINMENTS_EMPTY, CONTAINMENTS_SIMPLE, CONTAINMENTS_COMPLEX, \
    PATHS_EMPTY, PATHS_SIMPLE, PATHS_COMPLEX, \
    GFA1_EMPTY_FN, GFA1_SIMPLE_FN, GFA1_COMPLEX_FN, \
    GFA1_COMPLEX_SOFT_FN, GFA1_COMPLEX_HARD_FN

class TestReadGFA1(TestCase):
    """Tests for exfi.io.read_gfa.read_gfa1"""

    def test_empty(self):
        """exfi.io.read_gfa.read_gfa1: empty case"""
        gfa1 = read_gfa1(GFA1_EMPTY_FN)
        self.assertTrue(gfa1['header'].equals(HEADER))
        self.assertTrue(gfa1['segments'].equals(SEGMENTS_EMPTY))
        self.assertTrue(gfa1['links'].equals(LINKS_EMPTY))
        self.assertTrue(gfa1['containments'].equals(CONTAINMENTS_EMPTY))
        self.assertTrue(gfa1['paths'].equals(PATHS_EMPTY))

    def test_simple(self):
        """exfi.io.read_gfa.read_gfa1: simple case"""
        gfa1 = read_gfa1(GFA1_SIMPLE_FN)
        self.assertTrue(gfa1['header'].equals(HEADER))
        self.assertTrue(gfa1['segments'].equals(SEGMENTS_SIMPLE))
        self.assertTrue(gfa1['links'].equals(LINKS_SIMPLE))
        self.assertTrue(gfa1['containments'].equals(CONTAINMENTS_SIMPLE))
        self.assertTrue(gfa1['paths'].equals(PATHS_SIMPLE))

    def test_complex(self):
        """exfi.io.read_gfa.read_gfa1: complex case"""
        gfa1 = read_gfa1(GFA1_COMPLEX_FN)
        self.assertTrue(gfa1['header'].equals(HEADER))
        self.assertTrue(gfa1['segments'].equals(SEGMENTS_COMPLEX))
        self.assertTrue(gfa1['links'].equals(LINKS_COMPLEX))
        self.assertTrue(gfa1['containments'].equals(CONTAINMENTS_COMPLEX))
        self.assertTrue(gfa1['paths'].equals(PATHS_COMPLEX))

    def test_complex_soft(self):
        """exfi.io.read_gfa.read_gfa1: complex and soft masking case"""
        gfa1 = read_gfa1(GFA1_COMPLEX_SOFT_FN)
        self.assertTrue(gfa1['header'].equals(HEADER))
        self.assertTrue(gfa1['segments'].equals(SEGMENTS_COMPLEX_SOFT))
        self.assertTrue(gfa1['links'].equals(LINKS_COMPLEX))
        self.assertTrue(gfa1['containments'].equals(CONTAINMENTS_COMPLEX))
        self.assertTrue(gfa1['paths'].equals(PATHS_COMPLEX))

    def test_complex_hard(self):
        """exfi.io.read_gfa.read_gfa1: complex and hard masking case"""
        gfa1 = read_gfa1(GFA1_COMPLEX_HARD_FN)
        self.assertTrue(gfa1['header'].equals(HEADER))
        self.assertTrue(gfa1['segments'].equals(SEGMENTS_COMPLEX_HARD))
        self.assertTrue(gfa1['links'].equals(LINKS_COMPLEX))
        self.assertTrue(gfa1['containments'].equals(CONTAINMENTS_COMPLEX))
        self.assertTrue(gfa1['paths'].equals(PATHS_COMPLEX))



if __name__ == '__main__':
    main()
