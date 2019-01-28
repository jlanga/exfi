#!/usr/bin/env python3

"""tests.test_io.test_gff3_to_bed.py: tests for exfi.io.gff3_to_bed.py"""

from unittest import TestCase, main

from exfi.io.gff3_to_bed import gff3_to_bed3

from tests.io.gff3 import \
    GFF3_EMPTY_FN, GFF3_ENSEMBL_FN, GFF3_GMAP_FN

from tests.io.bed import \
    BED3_EMPTY, BED3_ENSEMBL, BED3_GMAP




class TestGFF3ToBED(TestCase):
    """Tests for exfi.io.gff3_to_bed.gff3_to_bed3"""

    def test_empty(self):
        """exfi.io.gff3_to_bed.gff3_to_bed3: empty case"""
        observed = gff3_to_bed3(GFF3_EMPTY_FN)
        self.assertTrue(observed.equals(BED3_EMPTY))

    def test_ensembl(self):
        """exfi.io.gff3_to_bed.gff3_to_bed3: ensembl case"""
        observed = gff3_to_bed3(GFF3_ENSEMBL_FN, mode="ensembl")
        print("Observed", observed.values.tolist(), observed.dtypes, sep='\n')
        self.assertTrue(observed.equals(BED3_ENSEMBL))

    def test_gmap(self):
        """exfi.io.gff3_to_bed.gff3_to_bed3: gmap case"""
        observed = gff3_to_bed3(GFF3_GMAP_FN, mode="gmap")
        print("Observed", observed.values.tolist(), observed.dtypes, sep='\n')
        self.assertTrue(observed.equals(BED3_GMAP))

    # def test_ncbi(self):
    #     """exfi.io.gff3_to_bed.gff3_to_bed3: ncbi case"""
    #     observed = gff3_to_bed3(GFF3_NCBI_FN, mode=ncbi)
    #     self.assertTrue(observed.equals(BED3_NCBI))




if __name__ == '__main__':
    main()
