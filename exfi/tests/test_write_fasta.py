#!/usr/bin/env python3
from unittest import TestCase
from exfi.read_fasta import read_fasta
from exfi.write_fasta import write_fasta

class TestWriteFasta(TestCase):
    devnull = "/dev/null"
    empty_file = "exfi/tests/files/empty_file"

    one_seq = "exfi/tests/files/one_sequence.fa"
    one_seq_bis = "exfi/tests/files/one_sequence_bis.fa"

    two_seqs = "exfi/tests/files/two_sequences.fa"
    two_seqs_bis = "exfi/tests/files/two_sequences_bis.fa"

    def test_empty_file(self):
        """Write empty file"""
        fasta_iterator = read_fasta(self.devnull)

        write_fasta(
            generator_seq= fasta_iterator,
            filename = self.empty_file
        )
        
        actual = [x for x in read_fasta(self.empty_file)]
        expected = [x for x in read_fasta(self.devnull)]
        self.assertEqual(actual, expected)

    def test_one_sequence(self):
        """Write one sequence"""
        fasta_iterator = read_fasta(self.one_seq)

        write_fasta(
            generator_seq = fasta_iterator,
            filename = self.one_seq_bis
        )

        actual = [x for x in read_fasta(self.one_seq)][0]
        expected = [x for x in read_fasta(self.one_seq_bis)][0]
        self.assertEqual(actual.id, expected.id)
        self.assertEqual(actual.seq, expected.seq)

    def test_two_sequences(self):
        """Write two sequences"""
        fasta_iterator = read_fasta(self.two_seqs)

        write_fasta(
            generator_seq= fasta_iterator,
            filename= self.two_seqs_bis
        )
        actual = [x for x in read_fasta(self.two_seqs)]
        expected = [x for x in read_fasta(self.two_seqs_bis)]
        self.assertEqual(
            sorted([x.id for x in actual]),
            sorted([x.id for x in expected])
        )
        self.assertEqual(
            sorted([str(x.seq) for x in actual]),
            sorted([str(x.seq) for x in expected])
        )

    @classmethod
    def tearDownClass(self):
        import os
        os.remove(self.empty_file)
        os.remove(self.one_seq_bis)
        os.remove(self.two_seqs_bis)