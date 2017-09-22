#!/usr/bin/env python3


from unittest import TestCase
from Bio import SeqIO
from exfi.reduce_exons import reduce_exons
from exfi.tests.auxiliary_functions import \
    CustomAssertions, \
    _fasta_to_list

class TestReduceExons(TestCase, CustomAssertions):

    @classmethod
    def test_empty_sequence(self):
        """reduce_exons.py: test empty fasta"""
        fasta_in = "exfi/tests/files/reduce_exons/empty_sequence_input.fa"
        fasta_res = "exfi/tests/files/reduce_exons/empty_sequence_output.fa"
        self.assertEqualListOfSeqrecords(
            list(reduce_exons(_fasta_to_list(fasta_in))),
            _fasta_to_list(fasta_res)
        )

    @classmethod
    def test_one_exon(self):
        """reduce_exons.py: test single exon"""
        fasta_in = "exfi/tests/files/reduce_exons/one_exon_input.fa"
        fasta_res = "exfi/tests/files/reduce_exons/one_exon_output.fa"
        self.assertEqualListOfSeqrecords(
            list(reduce_exons(_fasta_to_list(fasta_in))),
            _fasta_to_list(fasta_res)
        )

    @classmethod
    def test_same_exon(self):
        """reduce_exons.py: read the same exon twice"""
        fasta_in = "exfi/tests/files/reduce_exons/same_exon_input.fa"
        fasta_res = "exfi/tests/files/reduce_exons/same_exon_output.fa"
        self.assertEqualListOfSeqrecords(
            list(reduce_exons(_fasta_to_list(fasta_in))),
            _fasta_to_list(fasta_res)
        )

    @classmethod
    def test_different_exons(self):
        """reduce_exons.py: read the same exon twice"""
        fasta_in = "exfi/tests/files/reduce_exons/different_exons_input.fa"
        fasta_res = "exfi/tests/files/reduce_exons/different_exons_output.fa"
        self.assertEqualListOfSeqrecords(
            list(reduce_exons(_fasta_to_list(fasta_in))),
            _fasta_to_list(fasta_res)
        )
