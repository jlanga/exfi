#!/usr/bin/env python3


from unittest import TestCase
from exfi.reduce_exons import reduce_exons
from tests.auxiliary_functions import \
    CustomAssertions, \
    _fasta_to_list

class TestReduceExons(TestCase, CustomAssertions):

    def test_empty_sequence(self):
        """reduce_exons: empty"""
        actual = list(reduce_exons(_fasta_to_list(
            "tests/files/reduce_exons/empty_sequence_input.fa"
        )))
        expected = _fasta_to_list(
            "tests/files/reduce_exons/empty_sequence_output.fa"
        )
        self.assertEqualListOfSeqrecords(actual, expected)

    def test_one_exon(self):
        """reduce_exons: single exon"""
        actual = list(reduce_exons(_fasta_to_list(
            "tests/files/reduce_exons/one_exon_input.fa"
        )))
        expected = _fasta_to_list(
            "tests/files/reduce_exons/one_exon_output.fa"
        )
        self.assertEqualListOfSeqrecords(actual, expected)

    def test_same_exon(self):
        """reduce_exons: same exon twice"""
        actual = list(reduce_exons(_fasta_to_list(
            "tests/files/reduce_exons/same_exon_input.fa"
        )))
        expected = _fasta_to_list(
            "tests/files/reduce_exons/same_exon_output.fa"
        )
        self.assertEqualListOfSeqrecords(actual, expected)

    def test_different_exons(self):
        """reduce_exons: different exons"""
        actual = list(reduce_exons(_fasta_to_list(
            "tests/files/reduce_exons/different_exons_input.fa"
        )))
        expected = _fasta_to_list(
            "tests/files/reduce_exons/different_exons_output.fa"
        )
        self.assertEqualListOfSeqrecords(actual, expected)
