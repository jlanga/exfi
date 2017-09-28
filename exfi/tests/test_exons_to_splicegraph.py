#!/usr/bin/env python3

import unittest
from exfi.exons_to_splicegraph import \
    exons_to_df
import pandas as pd
from Bio import SeqIO


class TestExonsToDF(unittest.TestCase):

    def test_empty_index(self):
        """exons_to_splicegraph.py,  check if an empty exome generates an empty
        DataFrame"""
        self.assertTrue(
            exons_to_df({})\
            .equals(
                pd.DataFrame(
                    columns=['transcript_id', 'start', 'end', 'exon_id', 'score', 'strand']
                )
            )
        )

    def test_one_entry(self):
        """exons_to_splicegraph.py,  single exon file to DataFrame."""
        self.assertTrue(
            exons_to_df(
                SeqIO.index(
                    filename="exfi/tests/files/exons_to_splicegraph/single.fa",
                    format="fasta"
                )
            )\
            .equals(
                pd.DataFrame(
                    data=[["ENSDART00000161035.1", 0, 326, "EXON00000000001", 0, '+']],
                    columns=['transcript_id', 'start', 'end', 'exon_id', 'score', 'strand']
                )
            )

        )

    def test_multiple(self):
        self.assertTrue(
            exons_to_df(
                SeqIO.index(
                    filename="exfi/tests/files/exons_to_splicegraph/different_transcripts.fa",
                    format="fasta"
                )
            )\
            .equals(
                pd.DataFrame(
                    data=[
                        ["ENSDART00000161035.1", 397, 472, "EXON00000000002", 0, "+"],
                        ["ENSDART00000165342.1", 1176, 1324, "EXON00000000015", 0, "+"],
                        ["ENSDART00000161035.1", 0, 326, "EXON00000000001", 0, "+"],
                        ["ENSDART00000165342.1", 125, 304, "EXON00000000005", 0, "+"],
                        ["ENSDART00000165342.1", 746, 851, "EXON00000000010", 0, "+"],
                        ["ENSDART00000165342.1", 974, 1097, "EXON00000000013", 0, "+"],
                        ["ENSDART00000165342.1", 854, 886, "EXON00000000011", 0, "+"],
                        ["ENSDART00000165342.1", 1098, 1175, "EXON00000000014", 0, "+"],
                        ["ENSDART00000165342.1", 5, 127, "EXON00000000004", 0, "+"],
                        ["ENSDART00000165342.1", 645, 746, "EXON00000000009", 0, "+"],
                        ["ENSDART00000165342.1", 317, 460, "EXON00000000006", 0, "+"],
                        ["ENSDART00000165342.1", 591, 650, "EXON00000000008", 0, "+"],
                        ["ENSDART00000165342.1", 459, 592, "EXON00000000007",0, "+"],
                        ["ENSDART00000165342.1", 899, 953, "EXON00000000012",0, "+"],
                        ["ENSDART00000161035.1", 477, 523, "EXON00000000003",0, "+"],
                    ],
                    columns=['transcript_id', 'start', 'end', 'exon_id', 'score', 'strand']
                )\
                .sort_values(['transcript_id', 'start', 'end'])
            )
        )
