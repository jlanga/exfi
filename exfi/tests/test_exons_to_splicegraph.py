#!/usr/bin/env python3

import unittest
from exfi.exons_to_splicegraph import \
    exons_to_df, \
    exon_to_coordinates, \
    transcript_to_path

import pandas as pd
import numpy as np
from Bio import SeqIO


class TestExonsToDF(unittest.TestCase):

    def test_empty_index(self):
        """exons_to_splicegraph.py:  check if an empty exome generates an empty
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
        """exons_to_splicegraph.py:  single exon file to DataFrame"""
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
        """exons_to_splicegraph.py: multiple transcript - multiple exon file to DataFrame."""
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

class TestExonsToCoordinates(unittest.TestCase):

    def test_empty(self):
        """exons_to_splicegraph.py: Get coordinates of an empty file"""
        self.assertEqual(
            exon_to_coordinates({}),
            {}
        )

    def test_single(self):
        """exons_to_splicegraph.py: Get coordinates of a single exon"""
        self.assertEqual(
            exon_to_coordinates(
                SeqIO.index(
                    filename="exfi/tests/files/exons_to_splicegraph/single.fa",
                    format="fasta"
                )
            ),
            {"EXON00000000001": [("ENSDART00000161035.1", 0, 326)]}
        )

    def test_multiple(self):
        """exons_to_splicegraph.py: Get coordinates of a single exon"""
        self.assertEqual(
            exon_to_coordinates(
                SeqIO.index(
                    filename="exfi/tests/files/exons_to_splicegraph/different_transcripts.fa",
                    format="fasta"
                )
            ),
            {
                "EXON00000000002": [("ENSDART00000161035.1", 397, 472)],
                "EXON00000000015": [("ENSDART00000165342.1", 1176, 1324)],
                "EXON00000000001": [("ENSDART00000161035.1", 0, 326)],
                "EXON00000000005": [("ENSDART00000165342.1", 125, 304)],
                "EXON00000000010": [("ENSDART00000165342.1", 746, 851)],
                "EXON00000000013": [("ENSDART00000165342.1", 974, 1097)],
                "EXON00000000011": [("ENSDART00000165342.1", 854, 886)],
                "EXON00000000014": [("ENSDART00000165342.1", 1098, 1175)],
                "EXON00000000004": [("ENSDART00000165342.1", 5, 127)],
                "EXON00000000009": [("ENSDART00000165342.1", 645, 746)],
                "EXON00000000006": [("ENSDART00000165342.1", 317, 460)],
                "EXON00000000008": [("ENSDART00000165342.1", 591, 650)],
                "EXON00000000007": [("ENSDART00000165342.1", 459, 592)],
                "EXON00000000012": [("ENSDART00000165342.1", 899, 953)],
                "EXON00000000003": [("ENSDART00000161035.1", 477, 523)]
            }
        )


class TestTranscriptToPath(unittest.TestCase):

    def test_empty(self):
        """exons_to_splicegraph.py: convert an empty exome to path"""
        self.assertTrue(
            transcript_to_path(
                exons_to_df({})
            )\
            .equals(
                pd.DataFrame(
                    data= np.empty((0,2), dtype=np.float64),
                    columns=['transcript_id', 'path']
                )\
                .set_index('transcript_id')
            )
        )


    def test_single(self):
        """exons_to_splicegraph.py: convert an single exon transcript to path"""
        self.assertTrue(
            transcript_to_path(
                exons_to_df(
                    SeqIO.index(
                        filename="exfi/tests/files/exons_to_splicegraph/single.fa",
                        format="fasta"
                    )
                )
            )\
            .equals(
                pd.DataFrame(
                    data=[["ENSDART00000161035.1", ["EXON00000000001"]]],
                    columns=['transcript_id', 'path']
                )\
                .set_index('transcript_id')
            )
        )

    def test_multiple(self):
        """exons_to_splicegraph.py: convert an single exon transcript to path"""
        self.assertTrue(
            transcript_to_path(
                exons_to_df(
                    SeqIO.index(
                        filename="exfi/tests/files/exons_to_splicegraph/different_transcripts.fa",
                        format="fasta"
                    )
                )
            )\
            .equals(
                pd.DataFrame(
                    data=[
                        ["ENSDART00000161035.1",
                            ["EXON00000000001", "EXON00000000002", "EXON00000000003"]
                        ],
                        ["ENSDART00000165342.1",
                            ["EXON00000000004", "EXON00000000005", "EXON00000000006",
                            "EXON00000000007", "EXON00000000008", "EXON00000000009",
                            "EXON00000000010", "EXON00000000011", "EXON00000000012",
                            "EXON00000000013", "EXON00000000014", "EXON00000000015"]
                        ]
                    ],
                    columns=['transcript_id', 'path']
                )\
                .set_index('transcript_id')
            )
        )


if __name__ == '__main__':
    unittest.main()
