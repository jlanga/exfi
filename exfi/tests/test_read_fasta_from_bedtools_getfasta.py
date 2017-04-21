#!/usr/bin/env python3

from unittest import TestCase

from Bio import SeqIO
from exfi.read_fasta_from_bedtools_getfasta import read_fasta_from_bedtools_getfasta

class TestReadFastaFromBedtoolsGetfasta(TestCase):

    def test_empty_sequence(self):
        """Test empty"""
        file = "/dev/null"
        records = [x for x in list(SeqIO.parse(handle= file, format= "fasta"))]
        actual = [x for x in read_fasta_from_bedtools_getfasta(records)]
        expected = []
        self.assertEqual(actual, expected)

    def test_one_sequence(self):
        """Test one sequence"""
        file = "exfi/tests/files/one_sequence.fa"
        records = [x for x in list(SeqIO.parse(handle= file, format= "fasta"))]
        actual = [x.id for x in read_fasta_from_bedtools_getfasta(records)]
        expected = ["ENSDARE00000830915:12-19;1"]
        self.assertEqual(actual, expected)

    def test_two_exons_one_transcript(self):
        """Exons from same transcript"""
        file = "exfi/tests/files/two_seqs_same_transcript.fa"
        records = [x for x in list(SeqIO.parse(handle= file, format= "fasta"))]
        actual = [x.id for x in read_fasta_from_bedtools_getfasta(records)]
        expected = ["ENSDARE00000830915:12-31;1", "ENSDARE00000830915:17-25;1"]
        self.assertEqual(actual, expected)

    def test_one_exon_two_transcripts(self):
        """Exons from different transcripts"""
        file = "exfi/tests/files/two_seqs_different_transcript.fa"
        records = [x for x in list(SeqIO.parse(handle= file, format= "fasta"))]
        actual = [x.id for x in read_fasta_from_bedtools_getfasta(records)]
        expected = ['test1:15-25;1', 'test2:129-129;1']
        self.assertEqual(actual, expected)

    def test_multiple_exons_multiple_transcripts(self):
        """Test a realistic case"""
        file = "exfi/tests/files/multiple_exons_multiple_transcripts.fa"
        records = [x for x in list(SeqIO.parse(handle= file, format= "fasta"))]
        actual = [x.id for x in read_fasta_from_bedtools_getfasta(records)]
        expected = [
            'test1:15-29;1', 'test1:31-51;1', 'test1:60-70;1',
            'test2:15-30;1', 'test3:12-25;1', 'test3:24-37;1'
        ]
        self.assertEqual(actual, expected)
