#!/usr/bin/env python3


from unittest import TestCase
from exfi.extend import extend_left, extend_right
from Bio import SeqIO


class TestExtendLeft(TestCase):

    def test_extend_left_empty(self):
        """Process an empty list"""
        records = []
        actual = [extend_left(x, 27) for x in records]
        expected = []
        self.assertEqual(actual, expected)

    def test_extend_left_one(self):
        """Process just one sequence"""
        records = SeqIO.parse(
            handle="exfi/tests/files/one_sequence.fa",
            format="fasta"
        )
        actual = list(extend_left(records, 27))
        expected_ids = [
            "ENSDARE00000830915:12-19l" + letter for letter in "ACGT"
        ]
        expected_sequences = [  # seq was
            # GTAAGCCGCGGCGGTGTGTGTGTGTGTGTGTGTTCTCCGTCATCTGTGTTCTGCTGAATG
            letter + "GTAAGCCGCGGCGGTGTGTGTGTGTG" for letter in "ACGT"
        ]
        self.assertEqual([x.id for x in actual], expected_ids)
        self.assertEqual([str(x.seq) for x in actual], expected_sequences)

    def test_extend_left_multiple(self):
        """Process multiple sequences"""
        records = SeqIO.parse(
            handle="exfi/tests/files/four_sequences.fa",
            format="fasta"
        )
        actual = list(extend_left(records, 27))

        expected_ids = \
            ["tr1_e1l" + letter for letter in "ACGT"] + \
            ["tr1_e2l" + letter for letter in "ACGT"]
        self.assertEqual(
            [x.id for x in actual],
            expected_ids
        )

        expected_sequences = \
            [letter + "CATCGATCGATCAGTAGCTAGCTAGT" for letter in "ACGT"] + \
            [letter + "GCATGCTAGCTCATAGCTCAGTACGA" for letter in "ACGT"]
        self.assertEqual(
            [str(x.seq) for x in actual],
            expected_sequences
        )

        print([x.id for x in actual])
        print(expected_ids)


class TestExtendRight(TestCase):

    def test_extend_right_empty(self):
        """Process an empty list"""
        records = []
        actual = [extend_right(x, 27) for x in records]
        expected = []
        self.assertEqual(actual, expected)

    def test_extend_right_one(self):
        """Process just one sequence"""
        records = SeqIO.parse(
            handle="exfi/tests/files/one_sequence.fa",
            format="fasta"
        )
        actual = list(extend_right(records, 27))
        expected_ids = [
            "ENSDARE00000830915:12-19r" + letter for letter in "ACGT"
        ]
        expected_sequences = [  # seq was
            # GTAAGCCGCGGCGGTGTGTGTGTGTGTGTGTGTTCTCCGTCATCTGTGTTCTGCTGAATG
            "CTCCGTCATCTGTGTTCTGCTGAATG" + letter for letter in "ACGT"
        ]
        self.assertEqual([x.id for x in actual], expected_ids)
        self.assertEqual([str(x.seq) for x in actual], expected_sequences)

    def test_extend_right_multiple(self):
        """Process multiple sequences"""
        records = SeqIO.parse(
            handle="exfi/tests/files/four_sequences.fa",
            format="fasta"
        )
        actual = list(extend_right(records, 27))
        expected_ids = \
            ["tr1_e1r" + letter for letter in "ACGT"] + \
            ["tr1_e2r" + letter for letter in "ACGT"]
        self.assertEqual(
            [x.id for x in actual],
            expected_ids
        )

        expected_sequences = \
            ["CGATCGATCAGTAGCTAGCTAGTCAC" + letter for letter in "ACGT"] + \
            ["AGCTCAGTACGATCGACTACGAGTCA" + letter for letter in "ACGT"]

        print(str(x.seq) for x in actual)
        print()
        print([x for x in expected_sequences])
        self.assertEqual(
            [str(x.seq) for x in actual],
            expected_sequences
        )
