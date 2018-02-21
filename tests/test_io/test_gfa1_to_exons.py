#!/usr/bin/env python3

"""
Tests for exfi.io.gfa1_to_exons
"""

from unittest import TestCase, main

import filecmp
import tempfile
import os

from exfi.io.gfa1_to_exons import \
    gfa1_to_exons

from exfi.io.fasta_to_dict import \
    fasta_to_dict

from tests.data import \
    GFA_EMPTY_FN, GFA_SIMPLE_FN, GFA_COMPLEX_FN

EXONS_EMPTY_FN = "tests/io/exons_empty.fa"
EXONS_SIMPLE_FN = "tests/io/exons_simple.fa"
EXONS_COMPLEX_FN = "tests/io/exons_complex.fa"
EXONS_COMPLEX_SOFT_FN = "tests/io/exons_complex_soft.fa"
EXONS_COMPLEX_HARD_FN = "tests/io/exons_complex_hard.fa"

EXONS_EMPTY_DICT = fasta_to_dict(EXONS_EMPTY_FN)
EXONS_SIMPLE_DICT = fasta_to_dict(EXONS_SIMPLE_FN)
EXONS_COMPLEX_DICT = fasta_to_dict(EXONS_COMPLEX_FN)
EXONS_COMPLEX_SOFT_DICT = fasta_to_dict(EXONS_COMPLEX_SOFT_FN)
EXONS_COMPLEX_HARD_DICT = fasta_to_dict(EXONS_COMPLEX_HARD_FN)

class TestGFA1ToExons(TestCase):
    """Tests for gfa1_to_exons"""

    def test_empty(self):
        """exfi.io.gfa1_to_exons.gfa1_to_exons: empty case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(gfa_in_fn=GFA_EMPTY_FN, fasta_out_fn=tmp_file, masking="none")
        self.assertTrue(filecmp.cmp(tmp_file, EXONS_EMPTY_FN))
        os.remove(tmp_file)


    def test_simple(self):
        """exfi.io.gfa1_to_exons.gfa1_to_exons: simple case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(gfa_in_fn=GFA_SIMPLE_FN, fasta_out_fn=tmp_file)
        self.assertTrue(filecmp.cmp(tmp_file, EXONS_SIMPLE_FN))
        os.remove(tmp_file)

    def test_multiple(self):
        """exfi.io.gfa1_to_exons.gfa1_to_exons: complex case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(gfa_in_fn=GFA_COMPLEX_FN, fasta_out_fn=tmp_file)
        self.assertTrue(filecmp.cmp(tmp_file, EXONS_COMPLEX_FN))
        os.remove(tmp_file)

    def test_multiple_soft(self):
        """exfi.io.gfa1_to_exons.gfa1_to_exons: complex case and soft masking"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=GFA_COMPLEX_FN,
            fasta_out_fn=tmp_file,
            masking="soft"
        )
        self.assertTrue(
            filecmp.cmp(tmp_file, EXONS_COMPLEX_SOFT_FN)
        )
        os.remove(tmp_file)

    def test_multiple_hard(self):
        """exfi.io.gfa1_to_exons.gfa1_to_exons: complex case and hard masking"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=GFA_COMPLEX_FN,
            fasta_out_fn=tmp_file,
            masking="hard"
        )
        self.assertTrue(
            filecmp.cmp(tmp_file, EXONS_COMPLEX_HARD_FN)
        )
        os.remove(tmp_file)

if __name__ == '__main__':
    main()
