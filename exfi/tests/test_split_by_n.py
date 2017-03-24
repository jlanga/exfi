#!/usr/bin/env python3
from unittest import TestCase
from exfi.read_fasta import read_fasta
from exfi.split_by_n import split_by_n

class TestSplitByN(TestCase):
    """Test split_by_n function"""
    monoseq_file = "exfi/tests/files/one_sequence.fa"
    biseq_file = "exfi/tests/files/biseq.fa"
    multiseq_file = "exfi/tests/files/multiseq.fa"

    def test_split_empty_sequence(self):
        """Split an empty list"""
        iterator = read_fasta("/dev/null")
        actual = list(split_by_n(iterator))
        expected = []
        self.assertEqual(len(actual), len(expected))

    def test_split_monoexon_monoseq(self):
        """Split a list of only one sequence"""
        iterator = read_fasta(self.monoseq_file)
        actual = [x for x in split_by_n(iterator)]
        expected = [x for x in read_fasta(self.monoseq_file)]
        self.assertEqual(len(actual), len(expected))

    def test_split_multiexon_monoseq(self):
        """Split one sequence with Ns"""
        iterator = read_fasta(self.biseq_file)
        actual = len(list(split_by_n(iterator)))
        expected = 2
        self.assertEqual(actual, expected)

    def test_split_multiexon_multiseq(self):
        """Split a file of mutiple seqs with Ns"""
        iterator = read_fasta(self.multiseq_file)
        actual = len(list(split_by_n(iterator)))
        expected = 5
        self.assertEqual(actual, expected)