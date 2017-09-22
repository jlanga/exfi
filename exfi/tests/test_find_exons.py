#!/usr/bin/env python3

import unittest
from exfi.find_exons import \
    _process_output, \
    _get_fasta, \
    _find_exons_pipeline

from exfi.build_baited_bloom_filter import \
    _get_build_bf_command

from subprocess import Popen, PIPE
from Bio import SeqIO
from exfi.tests.auxiliary_functions import CustomAssertions
import tempfile
import shutil

def command_to_list(command):
    """Execute command and return output as list of strings"""
    process = Popen(command, stdout=PIPE, shell=False)
    results = list(_process_output(process))
    return results

class TestProcessOutput(unittest.TestCase):

    def test_empty_process(self):
        """find_exons.py: process an empty stream"""
        results = command_to_list(["cat", "/dev/null"])
        self.assertEqual(first=results, second=[])

    def test_simple_process(self):
        """find_exons.py: process an simple stream"""
        results = command_to_list(
            ["cat", "exfi/tests/files/find_exons/simple.bed"]
        )
        self.assertEqual(
            first=results,
            second=[("test", 4, 27)]
        )

    def test_big_process(self):
        """find_exons.py: process an big stream"""
        results = command_to_list(
            ["cat", "exfi/tests/files/find_exons/big.bed"]
        )
        self.assertEqual(
            first=results,
            second=[
                ("test1", 14, 27), ("test2", 15, 19),
                ("test2", 17, 21), ("test3", 19, 25)
            ]
        )


def fasta_to_dict(filename):
    """SeqIO.index wrapper for fasta files"""
    return SeqIO.index(filename=filename, format="fasta")


def fasta_to_list(filename):
    """SeqIO.parse wrapper for fasta files"""
    return list(SeqIO.parse(handle=filename, format="fasta"))

def getfasta_to_list(transcriptome_dict, iterable_of_bed):
    """Convert to a list the generator from getfasta"""
    return list(_get_fasta(transcriptome_dict, iterable_of_bed))


class TestGetFasta(unittest.TestCase, CustomAssertions):

    def test_empty_sequence_empty_bed(self):
        """find_exons.py: process an empty fasta and an empty bed"""
        transcriptome_dict={}
        iterable_of_bed=[]
        self.assertEqual(
            getfasta_to_list(transcriptome_dict, iterable_of_bed),
            []
        )

    def test_empty_sequence_one_bed(self):
        """find_exons.py: process an empty fasta and an empty bed"""
        transcriptome_dict={}
        iterable_of_bed=[("test1", 14, 27)]
        self.assertEqualListOfSeqrecords(
            getfasta_to_list(transcriptome_dict, iterable_of_bed),
            []
        )

    def test_one_sequence_empty_bed(self):
        """find_exons.py: process a simple fasta and an empty bed"""
        transcriptome_dict = fasta_to_dict(
            "exfi/tests/files/find_exons/single_sequence.fa"
        )
        iterable_of_bed = []
        self.assertEqualListOfSeqrecords(
            getfasta_to_list(transcriptome_dict, iterable_of_bed),
            []
        )

    def test_one_sequence_one_bed(self):
        """find_exons.py: process an single fasta and a single bed record"""
        transcriptome_dict = fasta_to_dict(
            "exfi/tests/files/find_exons/one_sequence_one_bed_input.fa"
        )
        iterable_of_bed = [("test1", 0, 60)]
        self.assertEqualListOfSeqrecords(
            getfasta_to_list(transcriptome_dict, iterable_of_bed),
            fasta_to_list(
                "exfi/tests/files/find_exons/one_sequence_one_bed_output.fa"
            )
        )

    def test_multiple_sequences_multiple_beds(self):
        """find_exons.py: process an multiline fasta and multple bed"""
        transcriptome_dict = fasta_to_dict(
            "exfi/tests/files/find_exons/multiple_sequences_multiple_beds_input.fa",
        )
        iterable_of_bed = [
            ("test1", 0, 60), ("test2", 0, 40), ("test3", 10, 20)
        ]
        self.assertEqualListOfSeqrecords(
            getfasta_to_list(transcriptome_dict, iterable_of_bed),
            fasta_to_list(
                "exfi/tests/files/find_exons/multiple_sequences_multiple_beds_output.fa",
            )
        )



def _silent_popen(command):
    """Create a Popen with no stderr and stdout"""
    return Popen(command,
        stdout=open("/dev/null", 'w'),
        stderr=open("/dev/null", 'w'),
        shell=False
    )

def _bf_and_process(reads_fns, transcriptome_fn):
    """Build the BF and process the reads"""
    tmp_dir = tempfile.mkdtemp()
    tmp_bf = tmp_dir + "transcriptome_noreads.bf"
    command = _get_build_bf_command("30", "100M", "1", "1", tmp_bf, reads_fns)
    process = _silent_popen(command)
    process.wait()
    results = _find_exons_pipeline(
        kmer=30,
        bloom_filter_fn=tmp_bf,
        transcriptome_fn=transcriptome_fn,
        max_fp_bases=5
    )
    shutil.rmtree(tmp_dir)
    return list(results)

class TestFindExonsPipeline(unittest.TestCase):

    def test_notranscriptome_noreads(self):
        """find_exons.py: Process an empty transcriptome and an empty BF"""
        reads_fns = ["/dev/null"]
        transcriptome_fn = "/dev/null"
        results = _bf_and_process(reads_fns, transcriptome_fn)
        self.assertEqual(results, [])

    def test_transcriptome_noreads(self):
        """find_exons.py: Process a small transcriptome and an empty BF"""
        reads_fns = ["/dev/null"]
        transcriptome_fn = 'exfi/tests/files/find_exons/small_transcriptome.fa'
        results = _bf_and_process(reads_fns, transcriptome_fn)
        self.assertEqual(results, [])

    def test_small_data(self):
        """find_exons.py: Process an empty transcriptome and a small BF"""
        reads_fns = [
            'exfi/tests/files/find_exons/reads_1.fq',
            'exfi/tests/files/find_exons/reads_2.fq'
        ]
        transcriptome_fn = 'exfi/tests/files/find_exons/small_transcriptome.fa'
        results = _bf_and_process(reads_fns, transcriptome_fn)
        self.assertEqual(results, [])


if __name__ == "__main__":
    unittest.main()
