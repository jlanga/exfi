#!/usr/bin/env python3

from unittest import TestCase

from exfi.read_fasta import read_fasta
from exfi.paste_sequences_with_ns import paste_sequences_with_ns

class TestPasteSequencesWithNS(TestCase):

    monoseq_file = "exfi/tests/files/one_sequence.fa"
    paste_two_file = "exfi/tests/files/paste_two.fa"
    four_sequences_file = "exfi/tests/files/four_sequences.fa"
    
    def test_paste_empty_sequence(self):
        """Paste nothing"""
        iterator = read_fasta("/dev/null")
        actual = list(paste_sequences_with_ns(iterator))
        expected = []
        self.assertEqual(len(actual), len(expected))

    def test_paste_one_sequence(self):
        """Paste from one sequence"""
        iterator = read_fasta(self.monoseq_file)
        actual = [x for x in paste_sequences_with_ns(iterator)]
        expected = [x for x in read_fasta(self.monoseq_file)]
        print(actual)
        self.assertEqual(len(actual), len(expected))

    def test_paste_one_block(self):
        """Join two sequences"""
        iterator = read_fasta(self.paste_two_file)
        actual = [x for x in paste_sequences_with_ns(iterator)]
        expected = 1
        self.assertEqual(len(actual), expected)

    def test_paste_two_regions(self):
        """Join four sequences into two regions"""
        iterator = read_fasta(self.four_sequences_file)
        actual = [x for x in paste_sequences_with_ns(iterator)]
        expected = 2
        self.assertEqual(len(actual), expected)



