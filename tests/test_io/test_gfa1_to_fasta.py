#!/usr/bin/python3

"""tests.test_io.test_gfa1_to_fasta.py: tests for exfi.io.gfa1_to_fasta.py"""

from unittest import TestCase, main

from tempfile import mkstemp
import os
import filecmp

from exfi.io.gfa1_to_fasta import gfa1_to_fasta

from tests.io.gfa1 import \
    GFA1_EMPTY_FN, GFA1_SIMPLE_FN, GFA1_COMPLEX_FN

from tests.io.fasta import \
    EXONS_EMPTY_FN, EXONS_SIMPLE_FN, EXONS_COMPLEX_FN, \
    EXONS_COMPLEX_SOFT_FN, EXONS_COMPLEX_HARD_FN, \
    GAPPED_EMPTY_FN, GAPPED_SIMPLE_FN, GAPPED_COMPLEX_FN, \
    GAPPED_COMPLEX_SOFT_FN, GAPPED_COMPLEX_HARD_FN


class TestGFA1ToExons(TestCase):
    """Tests for exfi.io.gfa1_to_fasta.gfa1_to_fasta, exon output"""

    def test_exons_empty(self):
        """exfi.io.gfa1_to_fasta.gfa1_to_fasta: empty exons case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        gfa1_to_fasta(fasta_out=tmp_file, gfa1_in=GFA1_EMPTY_FN)
        self.assertTrue(filecmp.cmp(tmp_file, EXONS_EMPTY_FN))
        os.remove(tmp_file)

    def test_exons_simple(self):
        """exfi.io.gfa1_to_fasta.gfa1_to_fasta: simple exon case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        gfa1_to_fasta(fasta_out=tmp_file, gfa1_in=GFA1_SIMPLE_FN)
        self.assertTrue(filecmp.cmp(tmp_file, EXONS_SIMPLE_FN))
        os.remove(tmp_file)

    def test_exons_complex(self):
        """exfi.io.gfa1_to_fasta.gfa1_to_fasta: complex exons case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        gfa1_to_fasta(fasta_out=tmp_file, gfa1_in=GFA1_COMPLEX_FN)
        self.assertTrue(filecmp.cmp(tmp_file, EXONS_COMPLEX_FN))
        os.remove(tmp_file)

    def test_exons_soft_masking(self):
        """exfi.io.gfa1_to_fasta.gfa1_to_fasta: complex and soft masking case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        gfa1_to_fasta(
            fasta_out=tmp_file,
            gfa1_in=GFA1_COMPLEX_FN,
            masking='soft'
        )
        self.assertTrue(filecmp.cmp(tmp_file, EXONS_COMPLEX_SOFT_FN))
        os.remove(tmp_file)

    def test_exons_hard_masking(self):
        """exfi.io.gfa1_to_fasta.gfa1_to_fasta: complex and hard masking case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        gfa1_to_fasta(
            fasta_out=tmp_file,
            gfa1_in=GFA1_COMPLEX_FN,
            masking='hard'
        )
        self.assertTrue(filecmp.cmp(tmp_file, EXONS_COMPLEX_HARD_FN))
        os.remove(tmp_file)



class TestGFA1ToGappedTranscripts(TestCase):
    """Tests for exfi.io.gfa1_to_fasta.gfa1_to_gapped_transcripts"""

    def test_transcripts_empty(self):
        """exfi.io.gfa1_to_fasta.gfa1_to_fasta: empty transcript case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        gfa1_to_fasta(
            fasta_out=tmp_file, gfa1_in=GFA1_EMPTY_FN, transcripts=True
        )
        self.assertTrue(filecmp.cmp(tmp_file, GAPPED_EMPTY_FN))
        os.remove(tmp_file)

    def test_simple(self):
        """exfi.io.gfa1_to_fasta.gfa1_to_fasta: simple transcript case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        gfa1_to_fasta(
            fasta_out=tmp_file, gfa1_in=GFA1_SIMPLE_FN, transcripts=True
        )
        self.assertTrue(filecmp.cmp(tmp_file, GAPPED_SIMPLE_FN))
        os.remove(tmp_file)

    def test_complex(self):
        """exfi.io.gfa1_to_fasta.gfa1_to_fasta: complex transcripts case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        gfa1_to_fasta(
            fasta_out=tmp_file, gfa1_in=GFA1_COMPLEX_FN, transcripts=True
        )
        self.assertTrue(filecmp.cmp(tmp_file, GAPPED_COMPLEX_FN))
        os.remove(tmp_file)

    def test_soft_masking(self):
        """exfi.io.gfa1_to_fasta.gfa1_to_gapped_transcripts: complex and soft masking case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        gfa1_to_fasta(
            fasta_out=tmp_file, gfa1_in=GFA1_COMPLEX_FN, masking='soft',
            transcripts=True
        )
        self.assertTrue(filecmp.cmp(tmp_file, GAPPED_COMPLEX_SOFT_FN))
        os.remove(tmp_file)

    def test_hard_masking(self):
        """exfi.io.gfa1_to_fasta.gfa1_to_gapped_transcripts: complex and hard masking case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        gfa1_to_fasta(
            fasta_out=tmp_file, gfa1_in=GFA1_COMPLEX_FN, masking='hard',
            transcripts=True
        )
        self.assertTrue(filecmp.cmp(tmp_file, GAPPED_COMPLEX_HARD_FN))
        os.remove(tmp_file)



if __name__ == '__main__':
    main()
