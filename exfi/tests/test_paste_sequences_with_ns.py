#!/usr/bin/env python3

from unittest import TestCase

from exfi.read_fasta import read_fasta
from exfi.paste_sequences_with_ns import paste_sequences_with_ns
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class TestPasteSequencesWithNS(TestCase):

    paste_one = "exfi/tests/files/paste_one.fa"
    paste_two_file = "exfi/tests/files/paste_two.fa"
    four_sequences_file = "exfi/tests/files/four_sequences.fa"
    
    def test_paste_empty_sequence(self):
        """Paste nothing"""
        iterator = read_fasta("/dev/null")
        actual = list(paste_sequences_with_ns(iterator))
        expected = []
        print(actual)
        print(expected)
        self.assertEqual(len(actual), len(expected))

    def test_paste_one_sequence(self):
        """Paste from one sequence"""
        iterator = read_fasta(self.paste_one)
        actual = [x for x in paste_sequences_with_ns(iterator)]
        expected = [SeqRecord(id="test1" , seq = Seq("ACCGTAGCATGCTAGCTACGTAGCTAGCTAGCTAG"))]
        print(actual)
        print(expected)
        self.assertEqual(
            [x.id for x in actual],
            [x.id for x in expected]
        )
        self.assertEqual(
            [x.seq for x in actual],
            [x.seq for x in expected]
        )

    def test_paste_one_block(self):
        """Join two sequences"""
        iterator = read_fasta(self.paste_two_file)
        actual = [x for x in paste_sequences_with_ns(iterator)]
        expected = 1
        print(actual)
        print(expected)
        self.assertEqual(len(actual), expected)

    def test_paste_two_regions(self):
        """Join four sequences into two regions"""
        iterator = read_fasta(self.four_sequences_file)
        actual = [x for x in paste_sequences_with_ns(iterator)]
        expected = 2
        print(actual)
        print(expected)
        self.assertEqual(len(actual), expected)
