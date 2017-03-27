#!/usr/bin/env python3

from unittest import TestCase

from exfi.read_fasta import read_fasta
from exfi.read_fasta_from_bedtools_getfasta import read_fasta_from_bedtools_getfasta

class TestReadFastaFromBedtoolsGetfasta(TestCase):

    def test_empty_sequence(self):
        """Test empty"""
        file = "/dev/null"
        records = [x for x in read_fasta(file)]
        actual = [x for x in read_fasta_from_bedtools_getfasta(records)]
        expected = []
        self.assertEqual(actual, expected)

    def test_one_sequence(self):
        """Test one sequence"""
        file = "exfi/tests/files/one_sequence.fa"
        records = [x for x in read_fasta(file)]
        actual = [x.id for x in read_fasta_from_bedtools_getfasta(records)]
        expected = ["ENSDARE00000830915_e1"]
        self.assertEqual(actual, expected)

    def test_two_exons_one_transcript(self):
        """Exons from same transcript"""
        file = "exfi/tests/files/two_seqs_same_transcript.fa"
        records = [x for x in read_fasta(file)]
        actual = [x.id for x in read_fasta_from_bedtools_getfasta(records)]
        expected = ["ENSDARE00000830915_e1", "ENSDARE00000830915_e2"]
        self.assertEqual(actual, expected)

    def test_one_exon_two_transcripts(self):
        """Exons from different transcripts"""
        file = "exfi/tests/files/two_seqs_different_transcript.fa"
        records = [x for x in read_fasta(file)]
        actual = [x.id for x in read_fasta_from_bedtools_getfasta(records)]
        expected = ["test1_e1", "test2_e1"]
        self.assertEqual(actual, expected)

    def test_multiple_exons_multiple_transcripts(self):
        """Test a realistic case"""
        file = "exfi/tests/files/multiple_exons_multiple_transcripts.fa"
        records = [x for x in read_fasta(file)]
        actual = [x.id for x in read_fasta_from_bedtools_getfasta(records)]
        expected = ["test1_e1", "test1_e2", "test1_e3", "test2_e1", "test3_e1", "test3_e2"]
        self.assertEqual(actual, expected)
