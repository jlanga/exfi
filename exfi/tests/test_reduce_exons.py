#!/usr/bin/env python3


from unittest import TestCase
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from exfi.reduce_exons import reduce_exons


class TestReduceExons(TestCase):

    def test_empty_sequence(self):
        """Test /dev/null"""
        records = SeqIO.parse(
            handle="/dev/null",
            format="fasta"
        )
        actual = list(reduce_exons(records))
        expected = []
        self.assertEqual(actual, expected)

    def test_one_exon(self):
        """Try a single exon"""
        records = list(SeqIO.parse(
            handle="exfi/tests/files/one_sequence.fa",
            format="fasta"
        ))
        actual = list(reduce_exons(records))
        expected = [SeqRecord(
            id="EXON00000000001",
            seq=Seq(
                "GTAAGCCGCGGCGGTGTGTGTGTGTGTGTGTGTTCTCCGTCATCTGTGTTCTGCTGAATG"
            )
        )]

        self.assertEqual(
            [record.id for record in actual],
            [record.id for record in expected]
        )

        self.assertEqual(
            [str(record.seq) for record in actual],
            [str(record.seq) for record in expected]
        )

    def test_same_exon(self):
        """Read the same exon twice"""
        records = list(SeqIO.parse(
            handle="exfi/tests/files/one_sequence.fa",
            format="fasta"
        )) * 2
        actual = list(reduce_exons(records))
        expected = [SeqRecord(
            id="EXON00000000001",
            seq=Seq(
                "GTAAGCCGCGGCGGTGTGTGTGTGTGTGTGTGTTCTCCGTCATCTGTGTTCTGCTGAATG"
            )
        )]
        self.assertEqual(
            [record.id for record in actual],
            [record.id for record in expected]
        )
        self.assertEqual(
            [str(record.seq) for record in actual],
            [str(record.seq) for record in expected]
        )

    def test_different_exons(self):
        """Read two different exons"""
        records = list(SeqIO.parse(
            handle="exfi/tests/files/two_sequences.fa",
            format="fasta"
        ))
        actual = list(reduce_exons(records))
        expected = [
            SeqRecord(
                id="EXON00000000001",
                seq=Seq(
                    "GTAAGCCGCGGCGGTGTGTGTGTGTGTGTGTGTTCTCCGTCAT"
                    "CTGTGTTCTGCTGAATG"
                )
            ),
            SeqRecord(
                id="EXON00000000002",
                seq=Seq(
                    "AGAATGCGAGTGAGTGTGTGCAGCCACAGTCGTCTGAGTTTCCT"
                    "GAAGGATTCTTCACGG"
                )
            )
        ]
        print(actual)
        print(expected)

        self.assertEqual(
            [record.id for record in actual],
            [record.id for record in expected]
        )
        self.assertEqual(
            [str(record.seq) for record in actual],
            [str(record.seq) for record in expected]
        )
