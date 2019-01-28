#!/usr/bin/env python3

"""tests.test_io.test_gfa1_to_bed.py: tests for exfi.io.gfa1_to_bed.py"""


from unittest import TestCase, main

from exfi.io.gfa1_to_bed import \
    gfa1_to_bed4

from tests.io.gfa1 import \
    GFA1_EMPTY_FN, GFA1_SIMPLE_FN, GFA1_COMPLEX_FN

from tests.io.bed import \
    BED4_EMPTY, BED4_SIMPLE, BED4_COMPLEX



class TestGFA1ToBED4(TestCase):
    """Tests for exfi.io.gfa1_to_bed.gfa1_to_bed4"""

    def test_empty(self):
        '''exfi.io.gfa1_to_bed.gfa1_to_bed4: empty case'''
        observed = gfa1_to_bed4(GFA1_EMPTY_FN)
        print("Observed", observed, observed.dtypes, sep="\n")
        print("Expected", BED4_EMPTY, BED4_EMPTY.dtypes, sep="\n")
        self.assertTrue(observed.equals(BED4_EMPTY))

    def test_simple(self):
        '''exfi.io.gfa1_to_bed.gfa1_to_bed4: simple case'''
        observed = gfa1_to_bed4(GFA1_SIMPLE_FN)
        self.assertTrue(observed.equals(BED4_SIMPLE))

    def test_complex(self):
        '''exfi.io.gfa1_to_bed.gfa1_to_bed4: complex case'''
        observed = gfa1_to_bed4(GFA1_COMPLEX_FN)
        self.assertTrue(observed.equals(BED4_COMPLEX))



if __name__ == '__main__':
    main()
