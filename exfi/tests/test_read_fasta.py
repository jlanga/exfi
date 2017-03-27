#!/usr/bin/env python3
from unittest import TestCase

from exfi.read_fasta import read_fasta

class TestReadFasta(TestCase):

    devnull = "/dev/null"
    one_seq = "exfi/tests/files/one_sequence.fa"
    two_seqs = "exfi/tests/files/two_sequences.fa"

    def test_empty_file(self):
        """Read an empty file"""
        _exhausted = object()
        actual = next(read_fasta(self.devnull), _exhausted)
        expected = _exhausted
        self.assertEqual(actual, expected)

    def test_one_sequence(self):
        """Read one sequence"""
        actual = [record for record in read_fasta(self.one_seq)][0]
        expected_id = "ENSDARE00000830915:12-19"
        expected_seq = "GTAAGCCGCGGCGGTGTGTGTGTGTGTGTGTGTTCTCCGTCATCTGTGTTCTGCTGAATG"
        self.assertEqual(actual.id, expected_id)
        self.assertEqual(str(actual.seq), expected_seq)

    def test_two_sequences(self):
        """Read two sequences"""
        actual = [record for record in read_fasta(self.two_seqs)]
        
        expected_id_1 = "ENSDARE00000830915"
        expected_seq_1 = "GTAAGCCGCGGCGGTGTGTGTGTGTGTGTGTGTTCTCCGTCATCTGTGTTCTGCTGAATG"

        expected_id_2 = "ENSDARE00000417255"
        expected_seq_2 = "AGAATGCGAGTGAGTGTGTGCAGCCACAGTCGTCTGAGTTTCCTGAAGGATTCTTCACGG"

        self.assertEqual(actual[0].id, expected_id_1)
        self.assertEqual(str(actual[0].seq), expected_seq_1)

        self.assertEqual(actual[1].id, expected_id_2)
        self.assertEqual(str(actual[1].seq), expected_seq_2)