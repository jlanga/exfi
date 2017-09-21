#!/usr/bin/env python3


from unittest import TestCase
from Bio import SeqIO
from exfi.reduce_exons import reduce_exons
from exfi.tests.auxiliary_functions import CustomAssertions


class TestReduceExons(TestCase, CustomAssertions):

    @classmethod
    def test_empty_sequence(self):
        """reduce_exons.py: test empty fasta"""
        records = SeqIO.parse(
            handle="exfi/tests/files/reduce_exons/empty_sequence_input.fa",
            format="fasta"
        )
        actual = reduce_exons(records)
        expected = SeqIO.parse(
            handle="exfi/tests/files/reduce_exons/empty_sequence_output.fa",
            format="fasta"
        )
        self.assertEqualListOfSeqrecords(
            list(actual),
            list(expected)
        )

    @classmethod
    def test_one_exon(self):
        """reduce_exons.py: test single exon"""
        records = SeqIO.parse(
            handle="exfi/tests/files/reduce_exons/one_exon_input.fa",
            format="fasta"
        )
        actual = reduce_exons(records)
        expected = SeqIO.parse(
            handle="exfi/tests/files/reduce_exons/one_exon_output.fa",
            format="fasta"
        )
        self.assertEqualListOfSeqrecords(
            list(actual),
            list(expected)
        )

    @classmethod
    def test_same_exon(self):
        """reduce_exons.py: read the same exon twice"""
        records = SeqIO.parse(
            handle="exfi/tests/files/reduce_exons/same_exon_input.fa",
            format="fasta"
        )
        actual = reduce_exons(records)
        expected = SeqIO.parse(
            handle="exfi/tests/files/reduce_exons/same_exon_output.fa",
            format="fasta"
        )
        self.assertEqualListOfSeqrecords(
            list(actual),
            list(expected)
        )

    @classmethod
    def test_different_exons(self):
        """reduce_exons.py: read the same exon twice"""
        records = SeqIO.parse(
            handle="exfi/tests/files/reduce_exons/different_exons_input.fa",
            format="fasta"
        )
        actual = reduce_exons(records)
        expected = SeqIO.parse(
            handle="exfi/tests/files/reduce_exons/different_exons_output.fa",
            format="fasta"
        )
        self.assertEqualListOfSeqrecords(
            list(actual),
            list(expected)
        )
