#!/usr/bin/env python3

from unittest import TestCase
from exfi.find_exons import \
    _process_output, \
    _get_fasta, \
    _find_exons_pipeline
    #find_exons
from subprocess import Popen, PIPE
from Bio import SeqIO
from exfi.tests.auxiliary_functions import CustomAssertions
import tempfile
import shutil

def command_to_list(command):
    """Execute command and return output as list of strings"""
    process = Popen(command, stdout=PIPE)
    results = list(_process_output(process))
    return results

class TestProcessOutput(TestCase):

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

class TestGetFasta(TestCase, CustomAssertions):

    def test_empty_sequence_empty_bed(self):
        """find_exons.py: process an empty fasta and an empty bed"""
        transcriptome_dict={}
        iterable_of_bed=[]
        results = list(_get_fasta(transcriptome_dict, iterable_of_bed))
        expected = []
        self.assertEqual(results, expected)

    def test_empty_sequence_one_bed(self):
        """find_exons.py: process an empty fasta and an empty bed"""
        transcriptome_dict={}
        iterable_of_bed=[("test1", 14, 27)]
        results = list(_get_fasta(transcriptome_dict, iterable_of_bed))
        expected = []
        self.assertEqual(results, expected)

    def test_one_sequence_empty_bed(self):
        """find_exons.py: process a simple fasta and an empty bed"""
        transcriptome_dict = fasta_to_dict(
            "exfi/tests/files/find_exons/single_sequence.fa"
        )
        iterable_of_bed = []
        results = list(_get_fasta(transcriptome_dict, iterable_of_bed))
        expected = []
        self.assertEqual(results, expected)

    def test_one_sequence_one_bed(self):
        """find_exons.py: process an single fasta and a single bed record"""
        transcriptome_dict = fasta_to_dict(
            "exfi/tests/files/find_exons/one_sequence_one_bed_input.fa"
        )
        iterable_of_bed = [("test1", 0, 60)]
        results = list(_get_fasta(transcriptome_dict, iterable_of_bed))
        expected = fasta_to_list(
            "exfi/tests/files/find_exons/one_sequence_one_bed_output.fa"
        )
        self.assertEqualListOfSeqrecords(results,expected)

    def test_multiple_sequences_multiple_beds(self):
        """find_exons.py: process an multiline fasta and multple bed"""
        transcriptome_dict = fasta_to_dict(
            "exfi/tests/files/find_exons/multiple_sequences_multiple_beds_input.fa",
        )
        iterable_of_bed = [
            ("test1", 0, 60), ("test2", 0, 40), ("test3", 10, 20)
        ]
        results = list(_get_fasta(transcriptome_dict, iterable_of_bed))
        expected = fasta_to_list(
            "exfi/tests/files/find_exons/multiple_sequences_multiple_beds_output.fa",
        )
        self.assertEqualListOfSeqrecords(results,expected)


class TestFindExonsPipeline(TestCase):

    def test_notranscriptome_noreads(self):
        """find_exons.py: Process an empty transcriptome and an empty BF"""
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/test.bf"
        process = Popen(
            ['abyss-bloom', 'build',
            '--kmer', "30",
            '--bloom-size', "100M",
            '--levels', "1",
            '--threads', "1",
            tmp_bf,
            '/dev/null'],
            stdout=open("/dev/null", 'w'),
            stderr=open("/dev/null", "w")
        )
        process.wait()
        results = _find_exons_pipeline(
            kmer=30,
            bloom_filter_fn=tmp_bf,
            transcriptome_fn="/dev/null",
            max_fp_bases=5
        )
        results = list(results)
        shutil.rmtree(tmp_dir)
        self.assertEqual(
            first=results,
            second=[]
        )

    def test_transcriptome_noreads(self):
        """find_exons.py: Process a small transcriptome and an empty BF"""
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/test.bf"
        process = Popen(['abyss-bloom', 'build',
                '--kmer', "30",
                '--bloom-size', "100M",
                '--levels', "1",
                '--threads', "1",
                tmp_bf,
                '/dev/null'],
            stdout=open('/dev/null', 'w'),
            stderr=open('/dev/null', 'w')
        )
        process.wait()
        results = _find_exons_pipeline(
            kmer=30,
            bloom_filter_fn=tmp_bf,
            transcriptome_fn='exfi/tests/files/find_exons/small_transcriptome.fa',
            max_fp_bases=5
        )
        results = list(results)
        shutil.rmtree(tmp_dir)
        self.assertEqual(
            first=results,
            second=[]
        )

    def test_small_data(self):
        """find_exons.py: Process an empty transcriptome and a small BF"""
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/test.bf"
        process = Popen(['abyss-bloom', 'build',
                '--kmer', "30",
                '--bloom-size', "100M",
                '--levels', "1",
                '--threads', "1",
                tmp_bf,
                'exfi/tests/files/find_exons/reads_1.fq',
                'exfi/tests/files/find_exons/reads_2.fq'],
            stdout=open('/dev/null', 'w'),
            stderr=open('/dev/null', 'w')
        )
        process.wait()
        results = _find_exons_pipeline(
            kmer=30,
            bloom_filter_fn=tmp_bf,
            transcriptome_fn='exfi/tests/files/find_exons/small_transcriptome.fa',
            max_fp_bases=5
        )
        results = list(results)
        shutil.rmtree(tmp_dir)
        self.assertEqual(
            first=results,
            second=[]
        )
