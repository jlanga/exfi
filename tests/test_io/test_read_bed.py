#!/usr/bin/env python3

"""tests.test_io.test_read_bed.py: tests for exfi.io.read_bed.py"""

from unittest import TestCase, main

from exfi.io.read_bed import read_bed3

from tests.io.bed import \
    BED3_EMPTY_FN, BED3_SIMPLE_FN, BED3_COMPLEX_FN, \
    BED3_EMPTY, BED3_SIMPLE, BED3_COMPLEX


class TestReadBed3(TestCase):
    """Tests for exfi.io.read_bed.read_bed3"""

    def test_empty(self):
        """exfi.io.read_bed.read_bed3: empty case"""
        observed = read_bed3(filename=BED3_EMPTY_FN)
        self.assertTrue(observed.equals(BED3_EMPTY))

    def test_simple(self):
        """exfi.io.read_bed.read_bed3: simple case"""
        observed = read_bed3(filename=BED3_SIMPLE_FN)
        self.assertTrue(observed.equals(BED3_SIMPLE))

    def test_complex(self):
        """exfi.io.read_bed.read_bed3: complex case"""
        observed = read_bed3(filename=BED3_COMPLEX_FN)
        self.assertTrue(observed.equals(BED3_COMPLEX))



if __name__ == '__main__':
    main()
