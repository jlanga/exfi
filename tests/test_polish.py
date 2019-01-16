#!/usr/bin/env python3

"""tests.test_polish_overlaps.py: Tests for the exfi.polish submodule"""

from unittest import \
    TestCase, \
    main

from exfi.polish import \
    polish_bed4

from tests.io.transcriptome_dicts import \
    TRANSCRIPTOME_EMPTY_DICT, \
    TRANSCRIPTOME_SIMPLE_DICT, \
    TRANSCRIPTOME_COMPLEX_DICT

from tests.io.bed import \
    BED4_EMPTY, BED4_SIMPLE, BED4_COMPLEX, \
    BED4_SIMPLE_POLISHED, BED4_COMPLEX_POLISHED



class TestPolishBED4(TestCase):
    """Tests for exfi.polish.polish_bed4"""


    def test_empty(self):
        """exfi.polish.polish_bed4: empty case"""
        observed = polish_bed4(BED4_EMPTY, TRANSCRIPTOME_EMPTY_DICT)
        self.assertTrue(observed.shape == (0, 4))

    def test_simple(self):
        """exfi.polish.polish_bed4: simple case"""
        observed = polish_bed4(BED4_SIMPLE, TRANSCRIPTOME_SIMPLE_DICT)
        print("Observed:\n", observed)
        print("Expected:\n", BED4_SIMPLE_POLISHED)
        self.assertTrue(observed.equals(BED4_SIMPLE_POLISHED))

    def test_complex(self):
        """exfi.polish.polish_bed4: complex case"""
        observed = polish_bed4(BED4_COMPLEX, TRANSCRIPTOME_COMPLEX_DICT)
        print("Observed:\n", observed)
        print("Expected:\n", BED4_COMPLEX_POLISHED)
        self.assertTrue(observed.equals(BED4_COMPLEX_POLISHED))



if __name__ == '__main__':
    main()
