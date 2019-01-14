#!/usr/bin/env python3

"""tests.test_polish_overlaps.py: Tests for the exfi.polish submodule"""

from unittest import \
    TestCase, \
    main

import pandas as pd

from exfi.polish import \
    polish_bed4

from tests.data import \
    TRANSCRIPTOME_EMPTY_DICT, \
    TRANSCRIPTOME_SIMPLE_DICT, \
    TRANSCRIPTOME_COMPLEX_DICT


# Test data
BED4_EMPTY = pd.DataFrame(
    data=None, columns=["chrom", "chromStart", "chromEnd", "name"]
)

BED4_SIMPLE = pd.DataFrame(
    data=[("ENSDART00000161035.1", 0, 326, "ENSDART00000161035.1:0-326")],
    columns=["chrom", "chromStart", "chromEnd", "name"]
)

BED4_COMPLEX = pd.DataFrame(
    data=[
        ["ENSDART00000161035.1", 397, 472, "ENSDART00000161035.1:397-472"],
        ["ENSDART00000161035.1", 0, 326, "ENSDART00000161035.1:0-326"],
        ["ENSDART00000161035.1", 477, 523, "ENSDART00000161035.1:477-523"],
        ["ENSDART00000165342.1", 1176, 1324, "ENSDART00000165342.1:1176-1324"],
        ["ENSDART00000165342.1", 125, 304, "ENSDART00000165342.1:125-304"],
        ["ENSDART00000165342.1", 746, 851, "ENSDART00000165342.1:746-851"],
        ["ENSDART00000165342.1", 974, 1097, "ENSDART00000165342.1:974-1097"],
        ["ENSDART00000165342.1", 854, 886, "ENSDART00000165342.1:854-886"],
        ["ENSDART00000165342.1", 1098, 1175, "ENSDART00000165342.1:1098-1175"],
        ["ENSDART00000165342.1", 5, 127, "ENSDART00000165342.1:5-127"],
        ["ENSDART00000165342.1", 645, 746, "ENSDART00000165342.1:645-746"],
        ["ENSDART00000165342.1", 317, 460, "ENSDART00000165342.1:317-460"],
        ["ENSDART00000165342.1", 591, 650, "ENSDART00000165342.1:591-650"],
        ["ENSDART00000165342.1", 459, 592, "ENSDART00000165342.1:459-592"],
        ["ENSDART00000165342.1", 899, 953, "ENSDART00000165342.1:899-953"]
    ],
    columns=["chrom", "chromStart", "chromEnd", "name"]
)


BED4_SIMPLE_POLISHED = BED4_SIMPLE
BED4_COMPLEX_POLISHED = pd.DataFrame(
    data=[
        ["ENSDART00000161035.1", 397, 472, "ENSDART00000161035.1:397-472"],
        ["ENSDART00000161035.1", 0, 326, "ENSDART00000161035.1:0-326"],
        ["ENSDART00000161035.1", 477, 523, "ENSDART00000161035.1:477-523"],
        ["ENSDART00000165342.1", 1176, 1324, "ENSDART00000165342.1:1176-1324"],
        ["ENSDART00000165342.1", 125, 304, "ENSDART00000165342.1:125-304"],
        ["ENSDART00000165342.1", 746, 851, "ENSDART00000165342.1:746-851"],
        ["ENSDART00000165342.1", 974, 1097, "ENSDART00000165342.1:974-1097"],
        ["ENSDART00000165342.1", 854, 886, "ENSDART00000165342.1:854-886"],
        ["ENSDART00000165342.1", 1098, 1175, "ENSDART00000165342.1:1098-1175"],
        ["ENSDART00000165342.1", 5, 127, "ENSDART00000165342.1:5-127"],
        ["ENSDART00000165342.1", 645, 746, "ENSDART00000165342.1:645-746"],
        ["ENSDART00000165342.1", 317, 460, "ENSDART00000165342.1:317-460"],
        ["ENSDART00000165342.1", 591, 650, "ENSDART00000165342.1:591-650"],
        ["ENSDART00000165342.1", 459, 592, "ENSDART00000165342.1:459-592"],
        ["ENSDART00000165342.1", 899, 953, "ENSDART00000165342.1:899-953"]
    ],
    columns=["chrom", "chromStart", "chromEnd", "name"]
)





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
