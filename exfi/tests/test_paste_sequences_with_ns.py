#!/usr/bin/env python3

from unittest import TestCase
from exfi.paste_sequences_with_ns import paste_sequences_with_ns
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

class TestPasteSequencesWithNS(TestCase):

    paste_one = "exfi/tests/files/paste_one.fa"
    paste_two_file = "exfi/tests/files/paste_two.fa"
    four_sequences_file = "exfi/tests/files/four_sequences.fa"
    
    def test_paste_empty_sequence(self):
        """Paste nothing"""
        iterator = list(SeqIO.parse(handle= "/dev/null", format= "fasta"))
        actual = list(paste_sequences_with_ns(iterator))
        expected = []
        print(actual)
        print(expected)
        self.assertEqual(len(actual), len(expected))

    def test_paste_one_sequence(self):
        """Paste from one sequence"""
        iterator = list(SeqIO.parse(handle= self.paste_one, format= "fasta"))
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
        iterator = list(SeqIO.parse(handle= self.paste_two_file, format= "fasta"))
        actual = [x for x in paste_sequences_with_ns(iterator)]
        expected = 1
        print(actual)
        print(expected)
        self.assertEqual(len(actual), expected)

    def test_paste_two_regions(self):
        """Join four sequences into two regions"""
        iterator = list(SeqIO.parse(handle= self.four_sequences_file, format= "fasta"))
        actual = [x for x in paste_sequences_with_ns(iterator)]
        expected = 2
        print(actual)
        print(expected)
        self.assertEqual(len(actual), expected)
