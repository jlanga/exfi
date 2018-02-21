#!/usr/bin/env python3

"""
Tests for exfi.io.gfa1_to_gapped_transcripts
"""


from unittest import TestCase, main

import filecmp
import tempfile
import os

from exfi.io.gfa1_to_gapped_transcripts import \
    gfa1_to_gapped_transcripts

from tests.data import \
    GFA_EMPTY_FN, GFA_SIMPLE_FN, GFA_COMPLEX_FN


GAPPED_EMPTY_FN = "tests/io/gapped_empty.fa"
GAPPED_SIMPLE_FN = "tests/io/gapped_simple.fa"
GAPPED_COMPLEX_FN = "tests/io/gapped_complex.fa"
GAPPED_COMPLEX_SOFT_FN = "tests/io/gapped_complex_soft.fa"
GAPPED_COMPLEX_HARD_FN = "tests/io/gapped_complex_hard.fa"

class TestGFA1ToGappedTranscripts(TestCase):
    """Tests for gfa1_to_gapped_transcript"""

    def test_empty(self):
        """exfi.io.gfa1_to_gapped_transcripts.gfa1_to_gapped_transcripts: empty case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcripts(gfa_in=GFA_EMPTY_FN, fasta_out=tmp_file)
        self.assertTrue(
            filecmp.cmp(tmp_file, GAPPED_EMPTY_FN)
        )
        os.remove(tmp_file)

    def test_simple(self):
        """exfi.io.gfa1_to_gapped_transcripts.gfa1_to_gapped_transcripts: simple case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcripts(gfa_in=GFA_SIMPLE_FN, fasta_out=tmp_file)
        self.assertTrue(
            filecmp.cmp(tmp_file, GAPPED_SIMPLE_FN)
        )
        os.remove(tmp_file)

    def test_multiple(self):
        """exfi.io.gfa1_to_gapped_transcripts.gfa1_to_gapped_transcripts: complex case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcripts(gfa_in=GFA_COMPLEX_FN, fasta_out=tmp_file)
        self.assertTrue(filecmp.cmp(
            tmp_file,
            GAPPED_COMPLEX_FN
        ))
        os.remove(tmp_file)

    def test_multiple_soft(self):
        """exfi.io.gfa1_to_gapped_transcripts.gfa1_to_gapped_transcripts: complex case and soft
        masking"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcripts(
            gfa_in=GFA_COMPLEX_FN,
            fasta_out=tmp_file,
            masking="soft"
        )
        self.assertTrue(
            filecmp.cmp(tmp_file, GAPPED_COMPLEX_SOFT_FN)
        )
        os.remove(tmp_file)


    def test_multiple_hard(self):
        """exfi.io.gfa1_to_gapped_transcripts.gfa1_to_gapped_transcripts: complex case and hard
        masking"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcripts(
            gfa_in=GFA_COMPLEX_FN,
            fasta_out=tmp_file,
            masking="hard"
        )
        self.assertTrue(filecmp.cmp(tmp_file, GAPPED_COMPLEX_HARD_FN))
        os.remove(tmp_file)

if __name__ == '__main__':
    main()
