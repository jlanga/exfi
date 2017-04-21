#!/usr/bin/env python3
from unittest import TestCase
from Bio import SeqIO
from exfi.split_by_n import split_by_n

class TestSplitByN(TestCase):
    """Test split_by_n function"""
    monoseq_file = "exfi/tests/files/one_sequence.fa"
    biseq_file = "exfi/tests/files/biseq.fa"
    multiseq_file = "exfi/tests/files/multiseq.fa"

    def test_split_empty_sequence(self):
        """Split an empty list"""
        iterator = list(SeqIO.parse(
            handle= "/dev/null",
            format= "fasta"
        ))
        actual = list(split_by_n(iterator))
        expected = []
        self.assertEqual(len(actual), len(expected))

    def test_split_monoexon_monoseq(self):
        """Split a list of only one sequence"""
        iterator = list(SeqIO.parse(
            handle= self.monoseq_file,
            format= "fasta"
        ))
        actual = [x for x in split_by_n(iterator)]
        expected = [x for x in list(SeqIO.parse(handle= self.monoseq_file, format= "fasta"))]
        self.assertEqual(len(actual), len(expected))

    def test_split_multiexon_monoseq(self):
        """Split one sequence with Ns"""
        iterator = list(SeqIO.parse(
            handle= self.biseq_file,
            format= "fasta"
        ))
        actual = len(list(split_by_n(iterator)))
        expected = 2
        self.assertEqual(actual, expected)

    def test_split_multiexon_multiseq(self):
        """Split a file of mutiple seqs with Ns"""
        iterator = list(SeqIO.parse(
            handle= self.multiseq_file,
            format= "fasta"
        ))
        actual = len(list(split_by_n(iterator)))
        expected = 5
        self.assertEqual(actual, expected)