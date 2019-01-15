#!/usr/bin/python3

"""tests.test_io.test_gfa1_to_fasta.py: tests for exfi.io.gfa1_to_fasta.py"""

from unittest import TestCase, main

from tempfile import mkstemp
import os
import filecmp

from exfi.io.gfa1_to_fasta import \
    gfa1_to_exons, gfa1_to_gapped_transcripts


from tests.io.gfa1 import \
    GFA1_EMPTY_FN, GFA1_SIMPLE_FN, GFA1_COMPLEX_FN

from tests.io.fasta import \
    EXONS_EMPTY_FN, EXONS_SIMPLE_FN, EXONS_COMPLEX_FN, \
    GAPPED_EMPTY_FN, GAPPED_SIMPLE_FN, GAPPED_COMPLEX_FN


class TestGFA1ToExons(TestCase):
    """Tests for exfi.io.gfa1_to_fasta.gfa1_to_exons"""

    def test_empty(self):
        """exfi.io.gfa1_to_fasta.gfa1_to_exons: empty case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        gfa1_to_exons(fasta_out=tmp_file, gfa1_in=GFA1_EMPTY_FN)
        self.assertTrue(filecmp.cmp(tmp_file, EXONS_EMPTY_FN))
        os.remove(tmp_file)

    def test_simple(self):
        """exfi.io.gfa1_to_fasta.gfa1_to_exons: simple case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        gfa1_to_exons(fasta_out=tmp_file, gfa1_in=GFA1_SIMPLE_FN)
        self.assertTrue(filecmp.cmp(tmp_file, EXONS_SIMPLE_FN))
        os.remove(tmp_file)

    def test_complex(self):
        """exfi.io.gfa1_to_fasta.gfa1_to_exons: complex case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        gfa1_to_exons(fasta_out=tmp_file, gfa1_in=GFA1_COMPLEX_FN)
        self.assertTrue(filecmp.cmp(tmp_file, EXONS_COMPLEX_FN))
        os.remove(tmp_file)

    # def test_soft_masking(self):
    #     """exfi.io.gfa1_to_fasta.gfa1_to_exons: complex and soft masking case"""
    #     pass
    #
    # def test_complex(self):
    #     """exfi.io.gfa1_to_fasta.gfa1_to_exons: complex and hard masking case"""
    #     pass



class TestGFA1ToGappedTranscripts(TestCase):
    """Tests for exfi.io.gfa1_to_fasta.gfa1_to_gapped_transcripts"""

    def test_empty(self):
        """exfi.io.gfa1_to_fasta.gfa1_to_gapped_transcripts: empty case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        gfa1_to_gapped_transcripts(fasta_out=tmp_file, gfa1_in=GFA1_EMPTY_FN)
        self.assertTrue(filecmp.cmp(tmp_file, GAPPED_EMPTY_FN))
        os.remove(tmp_file)

    def test_simple(self):
        """exfi.io.gfa1_to_fasta.gfa1_to_gapped_transcripts: simple case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        gfa1_to_gapped_transcripts(fasta_out=tmp_file, gfa1_in=GFA1_SIMPLE_FN)
        self.assertTrue(filecmp.cmp(tmp_file, GAPPED_SIMPLE_FN))
        os.remove(tmp_file)

    def test_complex(self):
        """exfi.io.gfa1_to_fasta.gfa1_to_gapped_transcripts: complex case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        gfa1_to_gapped_transcripts(fasta_out=tmp_file, gfa1_in=GFA1_COMPLEX_FN)
        self.assertTrue(filecmp.cmp(tmp_file, GAPPED_COMPLEX_FN))
        os.remove(tmp_file)

    # def test_soft_masking(self):
    #     """exfi.io.gfa1_to_fasta.gfa1_to_gapped_transcripts: complex and soft masking case"""
    #     pass
    #
    # def test_complex(self):
    #     """exfi.io.gfa1_to_fasta.gfa1_to_gapped_transcripts: complex and hard masking case"""
    #     pass



if __name__ == '__main__':
    main()
