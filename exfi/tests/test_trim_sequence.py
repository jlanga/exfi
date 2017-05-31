#!/usr/bin/env python3


from unittest import TestCase
from Bio import SeqIO
from exfi.trim_sequence import trim_sequence


class TestTrimSequence(TestCase):

    def test_trim_empty(self):
        """Trim an empty iterable"""
        iterator = []
        actual = [x for x in trim_sequence(iterator, 0, 0)]
        expected = []
        self.assertEqual(actual, expected)

    def test_trim_too_short(self):
        """Trim an impossible number of bases"""
        iterator = list(SeqIO.parse(
            handle="exfi/tests/files/one_sequence.fa",
            format="fasta"
        ))
        actual = [x for x in trim_sequence(iterator, 100, 100)]
        expected = []
        self.assertEqual(actual, expected)

    def test_trim_left_exact_length(self):
        """Trim the length of the sequence on the right"""
        iterator = list(SeqIO.parse(
            handle="exfi/tests/files/one_sequence.fa",
            format="fasta"
        ))
        actual = [x for x in trim_sequence(iterator, 60, 0)]
        expected = []
        self.assertEqual(actual, expected)

    def test_trim_right_exact_length(self):
        """Trim the length of the sequence"""
        iterator = list(SeqIO.parse(
            handle="exfi/tests/files/one_sequence.fa",
            format="fasta"
        ))
        actual = [x for x in trim_sequence(iterator, 0, 60)]
        expected = []
        self.assertEqual(actual, expected)

    def test_trim_both_exact_length(self):
        """Trim the length of the sequence"""
        iterator = list(SeqIO.parse(
            handle="exfi/tests/files/one_sequence.fa",
            format="fasta"
        ))
        actual = [x for x in trim_sequence(iterator, 30, 30)]
        expected = []
        self.assertEqual(actual, expected)

    def test_trim_left(self):
        """Trim bases only on the left side"""
        iterator = list(SeqIO.parse(
            handle="exfi/tests/files/one_sequence.fa",
            format="fasta"
        ))
        actual = [x for x in trim_sequence(iterator, 3, 0)]
        expected = "AGCCGCGGCGGTGTGTGTGTGTGTGTGTGTTCTCCGTCATCTGTGTTCTGCTGAATG"
        self.assertEqual(str(actual[0].seq), expected)

    def test_trim_right(self):
        """Trim bases only on the right side"""
        iterator = list(SeqIO.parse(
            handle="exfi/tests/files/one_sequence.fa",
            format="fasta"
        ))
        actual = [x for x in trim_sequence(iterator, 0, 10)]
        expected = "GTAAGCCGCGGCGGTGTGTGTGTGTGTGTGTGTTCTCCGTCATCTGTGTT"
        self.assertEqual(str(actual[0].seq), expected)

    def test_trim_both_sides(self):
        """Trim from both sides in one sequence"""
        iterator = list(SeqIO.parse(
            handle="exfi/tests/files/one_sequence.fa",
            format="fasta"
        ))
        actual = [x for x in trim_sequence(iterator, 3, 10)]
        expected = "AGCCGCGGCGGTGTGTGTGTGTGTGTGTGTTCTCCGTCATCTGTGTT"
        self.assertEqual(str(actual[0].seq), expected)

    def test_trim_multiple(self):
        """Trim multiple sequences, leaving one out"""
        iterator = list(SeqIO.parse(
            handle="exfi/tests/files/multiseq.fa",
            format="fasta"
        ))
        actual = [x for x in trim_sequence(iterator, 4, 4)]
        expected = [
            "GCCGCGGCGGTGTGTGNNNNNNNNNNNCCGTCATCTGTGTTCTGCTG",
            "TGCGAGTGAGTGTGTGCAGCNNNNNNNNNAGTTTCCTGAAGGATTCTTC"
        ]
        print(actual)
        print(expected)
        self.assertEqual([str(x.seq) for x in actual], expected)
