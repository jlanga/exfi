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

from tests.test_data import \
    GFA_EMPTY_FN, GFA_SIMPLE_FN, GFA_COMPLEX_FN, \
    GAPPED_EMPTY_FN, GAPPED_SIMPLE_FN, \
    GAPPED_COMPLEX_FN, GAPPED_COMPLEX_HARD_FN, GAPPED_COMPLEX_SOFT_FN


class TestGFA1ToGappedTranscripts(TestCase):
    """Tests for gfa1_to_gapped_transcript"""

    def test_empty(self):
        """gfa1_to_gapped_transcripts: empty case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcripts(gfa_in=GFA_EMPTY_FN, fasta_out=tmp_file)
        self.assertTrue(
            filecmp.cmp(tmp_file, GAPPED_EMPTY_FN)
        )
        os.remove(tmp_file)

    def test_simple(self):
        """gfa1_to_gapped_transcripts: simple case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcripts(gfa_in=GFA_SIMPLE_FN, fasta_out=tmp_file)
        self.assertTrue(
            filecmp.cmp(tmp_file, GAPPED_SIMPLE_FN)
        )
        os.remove(tmp_file)

    def test_multiple(self):
        """gfa1_to_gapped_transcripts: complex case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcripts(gfa_in=GFA_COMPLEX_FN, fasta_out=tmp_file)
        self.assertTrue(filecmp.cmp(
            tmp_file,
            GAPPED_COMPLEX_FN
        ))
        os.remove(tmp_file)

    def test_multiple_soft(self):
        """gfa1_to_gapped_transcripts: complex case and soft masking"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcripts(
            gfa_in=GFA_COMPLEX_FN,
            fasta_out=tmp_file,
            soft_mask_overlaps=True)
        self.assertTrue(
            filecmp.cmp(tmp_file, GAPPED_COMPLEX_SOFT_FN)
        )
        os.remove(tmp_file)


    def test_multiple_hard(self):
        """gfa1_to_gapped_transcripts: complex case and hard masking"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcripts(
            gfa_in=GFA_COMPLEX_FN,
            fasta_out=tmp_file,
            hard_mask_overlaps=True
        )
        self.assertTrue(
            filecmp.cmp(tmp_file, GAPPED_COMPLEX_HARD_FN)
        )
        os.remove(tmp_file)

if __name__ == '__main__':
    main()
