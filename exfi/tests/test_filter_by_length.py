#!/usr/bin/env python3
from unittest import TestCase

from exfi.read_fasta import read_fasta
from exfi.filter_by_length import filter_by_length

class TestFilterByLength(TestCase):
    """Test the filter_by_length function"""
    monoseq_file = "exfi/tests/files/one_sequence.fa"
    biseq_file = "exfi/tests/files/two_sequences.fa"
    multiseq_file = "exfi/tests/files/multiseq.fa"

    def test_empty_sequence(self):
        """Test an empty sequence"""
        actual = [filter_by_length(record, 20) for record in read_fasta("/dev/null")]
        expected = []
        self.assertEqual(actual, expected)

    def test_single_short_sequence(self):
        """Test a short sequence"""
        records = read_fasta(self.monoseq_file)
        actual = list(filter_by_length(records, 100))
        print(actual)
        expected = []
        self.assertEqual(actual, expected)

    def test_single_long_sequence(self):
        """Test a long sequence"""
        records = read_fasta(self.monoseq_file)
        actual = list(filter_by_length(records, 10))
        expected = list(read_fasta(self.monoseq_file))
        self.assertEqual(actual[0].id, expected[0].id)
        self.assertEqual(actual[0].seq, expected[0].seq)

    def test_two_sequences(self):
        """Test a file containing two sequences"""
        records = list(read_fasta(self.biseq_file))
        actual = list(filter_by_length(records, 10))
        expected = records
        self.assertEqual(len(actual), len(expected))

    def test_multiple_sequences(self):
        """Test a fasta with multiple sequences"""
        records = list(read_fasta(self.multiseq_file))
        actual = list(filter_by_length(records, 10))
        self.assertEqual(len(actual), 2)
