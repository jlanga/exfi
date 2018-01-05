#!/usr/bin/env python3

import unittest
from tests.auxiliary_functions import \
    CustomAssertions, \
    _command_to_list, \
    _fasta_to_dict, \
    _fasta_to_list, \
    _getfasta_to_list, \
    _bf_and_process


class TestProcessOutput(unittest.TestCase):

    def test_empty_process(self):
        """_command_to_list: process an empty stream"""
        results = _command_to_list(["cat", "/dev/null"])
        self.assertEqual(first=results, second=[])

    def test_simple_process(self):
        """_command_to_list: process an simple stream"""
        results = _command_to_list(
            ["cat", "tests/files/find_exons/simple.bed"]
        )
        self.assertEqual(
            first=results,
            second=[("test", 4, 27)]
        )

    def test_big_process(self):
        """_command_to_list: process an big stream"""
        results = _command_to_list(
            ["cat", "tests/files/find_exons/big.bed"]
        )
        self.assertEqual(
            first=results,
            second=[
                ("test1", 14, 27), ("test2", 15, 19),
                ("test2", 17, 21), ("test3", 19, 25)
            ]
        )


class TestGetFasta(unittest.TestCase, CustomAssertions):

    def test_empty_sequence_empty_bed(self):
        """_getfasta_to_list: process an empty fasta and an empty bed"""
        transcriptome_dict={}
        iterable_of_bed=[]
        self.assertEqual(
            _getfasta_to_list(transcriptome_dict, iterable_of_bed),
            []
        )

    def test_empty_sequence_one_bed(self):
        """_getfasta_to_list: process an empty fasta and an empty bed"""
        transcriptome_dict={}
        iterable_of_bed=[("test1", 14, 27)]
        self.assertEqualListOfSeqrecords(
            _getfasta_to_list(transcriptome_dict, iterable_of_bed),
            []
        )

    def test_one_sequence_empty_bed(self):
        """_getfasta_to_list: process a simple fasta and an empty bed"""
        transcriptome_dict = _fasta_to_dict(
            "tests/files/find_exons/single_sequence.fa"
        )
        iterable_of_bed = []
        self.assertEqualListOfSeqrecords(
            _getfasta_to_list(transcriptome_dict, iterable_of_bed),
            []
        )

    def test_one_sequence_one_bed(self):
        """_getfasta_to_list: process an single fasta and a single bed record"""
        transcriptome_dict = _fasta_to_dict(
            "tests/files/find_exons/one_sequence_one_bed_input.fa"
        )
        iterable_of_bed = [("test1", 0, 60)]
        self.assertEqualListOfSeqrecords(
            _getfasta_to_list(transcriptome_dict, iterable_of_bed),
            _fasta_to_list(
                "tests/files/find_exons/one_sequence_one_bed_output.fa"
            )
        )

    def test_multiple_sequences_multiple_beds(self):
        """_getfasta_to_list: process an multiline fasta and multple bed"""
        transcriptome_dict = _fasta_to_dict(
            "tests/files/find_exons/multiple_sequences_multiple_beds_input.fa",
        )
        iterable_of_bed = [
            ("test1", 0, 60), ("test2", 0, 40), ("test3", 10, 20)
        ]
        self.assertEqualListOfSeqrecords(
            _getfasta_to_list(transcriptome_dict, iterable_of_bed),
            _fasta_to_list(
                "tests/files/find_exons/multiple_sequences_multiple_beds_output.fa",
            )
        )


class TestFindExonsPipeline(unittest.TestCase):

    def test_notranscriptome_noreads(self):
        """_bf_and_process: Process an empty transcriptome and an empty BF"""
        reads_fns = ["/dev/null"]
        transcriptome_fn = "/dev/null"
        results = _bf_and_process(reads_fns, transcriptome_fn)
        self.assertEqual(results, [])

    def test_transcriptome_noreads(self):
        """_bf_and_process: Process a small transcriptome and an empty BF"""
        reads_fns = ["/dev/null"]
        transcriptome_fn = 'tests/files/find_exons/small_transcriptome.fa'
        results = _bf_and_process(reads_fns, transcriptome_fn)
        self.assertEqual(results, [])

    def test_small_data(self):
        """_bf_and_process: Process an empty transcriptome and a small BF"""
        reads_fns = [
            'tests/files/find_exons/reads_1.fq',
            'tests/files/find_exons/reads_2.fq'
        ]
        transcriptome_fn = 'tests/files/find_exons/small_transcriptome.fa'
        results = _bf_and_process(reads_fns, transcriptome_fn)
        self.assertEqual(results, [])


if __name__ == "__main__":
    unittest.main()
