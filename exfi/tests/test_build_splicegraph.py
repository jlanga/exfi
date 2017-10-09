#!/usr/bin/env python3

import unittest
from exfi.build_splicegraph import \
    exons_to_df, \
    exon_to_coordinates, \
    transcript_to_path, \
    compute_edge_overlaps, \
    build_splicegraph

import networkx as nx
import pandas as pd
from Bio import SeqIO

from exfi.io import \
    _clean_index


index_simple = _clean_index(SeqIO.index(
    filename="exfi/tests/files/build_splicegraph/single.fa",
    format="fasta"
))

index_different = _clean_index(SeqIO.index(
    filename="exfi/tests/files/build_splicegraph/different_transcripts.fa",
    format="fasta"
))

path_simple = {"ENSDART00000161035.1": ["EXON00000000001"]}

path_different = {
    "ENSDART00000161035.1":
        ["EXON00000000001", "EXON00000000002", "EXON00000000003"]
    ,
    "ENSDART00000165342.1":
        ["EXON00000000004", "EXON00000000005", "EXON00000000006",
        "EXON00000000007", "EXON00000000008", "EXON00000000009",
        "EXON00000000010", "EXON00000000011", "EXON00000000012",
        "EXON00000000013", "EXON00000000014", "EXON00000000015"]
}

transcriptome_simple = _clean_index(SeqIO.index(
    filename="exfi/tests/files/build_splicegraph/transcriptome_simple.fa",
    format="fasta"
))

transcriptome_different = _clean_index(SeqIO.index(
    filename="exfi/tests/files/build_splicegraph/transcriptome_different.fa",
    format="fasta"
))

def _prepare_overlaps(exons):
    """Compute splicegraph prior the computation of overlaps"""
    splice_graph = nx.DiGraph()
    exon_df = exons_to_df(exons)
    exon2coord = exon_to_coordinates(exons)
    splice_graph.add_nodes_from(exon2coord.keys())
    nx.set_node_attributes(
        G=splice_graph,
        name='coordinates',
        values = exon2coord
    )
    transcript2path = transcript_to_path(exon_df)
    for path in transcript2path.values():
        splice_graph.add_path(path)
    return splice_graph




class TestExonsToDF(unittest.TestCase):

    def test_empty_index(self):
        """build_splicegraph.py:  check if an empty exome generates an empty
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
        """build_splicegraph.py:  single exon file to DataFrame"""
        self.assertTrue(
            exons_to_df(
                index_simple
            )\
            .equals(
                pd.DataFrame(
                    data=[["ENSDART00000161035.1", 0, 326, "EXON00000000001", 0, '+']],
                    columns=['transcript_id', 'start', 'end', 'exon_id', 'score', 'strand']
                )
            )
        )

    def test_multiple(self):
        """build_splicegraph.py: multiple transcript - multiple exon file to DataFrame."""
        self.assertTrue(
            exons_to_df(index_different)\
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
        """build_splicegraph.py: Get coordinates of an empty file"""
        self.assertEqual(
            exon_to_coordinates({}),
            {}
        )

    def test_single(self):
        """build_splicegraph.py: Get coordinates of a single exon"""
        self.assertEqual(
            exon_to_coordinates(
                index_simple
            ),
            {"EXON00000000001": [("ENSDART00000161035.1", 0, 326)]}
        )

    def test_multiple(self):
        """build_splicegraph.py: Get coordinates of a single exon"""
        self.assertEqual(
            exon_to_coordinates(index_different),
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
        """build_splicegraph.py: convert an empty exome to path"""
        self.assertEqual(
            transcript_to_path(exons_to_df({})),
            {}
        )


    def test_single(self):
        """build_splicegraph.py: convert an single exon transcript to path"""
        self.assertEqual(
            transcript_to_path(exons_to_df(index_simple)),
            path_simple
        )

    def test_multiple(self):
        """build_splicegraph.py: convert an single exon transcript to path"""
        self.assertEqual(
            transcript_to_path(exons_to_df(index_different)),
            path_different
        )


class TestComputeOverlaps(unittest.TestCase):

    def test_empty_exome(self):
        """build_splicegraph.py: compute the overlaps of an empty exome"""
        splice_graph = _prepare_overlaps({})
        overlaps = compute_edge_overlaps(splice_graph)
        self.assertEqual(overlaps, {})

    def test_single_exon(self):
        """build_splicegraph.py: compute the overlaps of a single exon exome"""
        splice_graph =  _prepare_overlaps(index_simple)
        overlaps = compute_edge_overlaps(splice_graph)
        self.assertEqual(overlaps, {})

    def test_multiple_exons(self):
        """build_splicegraph.py: compute the overlaps of a simple exome"""
        splice_graph =  _prepare_overlaps(index_different)
        overlaps = compute_edge_overlaps(splice_graph)
        self.assertEqual(
            overlaps,
            {
                ('EXON00000000001', 'EXON00000000002'): -71,
                ('EXON00000000002', 'EXON00000000003'): -5,
                ('EXON00000000004', 'EXON00000000005'): 2,
                ('EXON00000000005', 'EXON00000000006'): -13,
                ('EXON00000000006', 'EXON00000000007'): 1,
                ('EXON00000000007', 'EXON00000000008'): 1,
                ('EXON00000000008', 'EXON00000000009'): 5,
                ('EXON00000000009', 'EXON00000000010'): 0,
                ('EXON00000000010', 'EXON00000000011'): -3,
                ('EXON00000000011', 'EXON00000000012'): -13,
                ('EXON00000000012', 'EXON00000000013'): -21,
                ('EXON00000000013', 'EXON00000000014'): -1,
                ('EXON00000000014', 'EXON00000000015'): -1
            }
        )


class TestBuildSplicegraph(unittest.TestCase):

    def test_empty(self):
        exons = {}
        self.assertTrue(
            nx.is_isomorphic(
                build_splicegraph(exons),
                nx.DiGraph()
            )
        )

    def test_simple(self):
        exons = index_simple
        expected = nx.DiGraph()
        expected.add_nodes_from(['EXON00000000001'])
        self.assertTrue(
            nx.is_isomorphic(
                build_splicegraph(exons),
                expected
            )
        )

    def test_multiple(self):
        exons = index_different
        expected = nx.DiGraph()
        expected.add_nodes_from(
            [
                "EXON00000000001", "EXON00000000002", "EXON00000000003"
            ] + \
            [
                "EXON00000000004", "EXON00000000005", "EXON00000000006",
                "EXON00000000007", "EXON00000000008", "EXON00000000009",
                "EXON00000000010", "EXON00000000011", "EXON00000000012",
                "EXON00000000013", "EXON00000000014", "EXON00000000015"
            ]
        )
        #paths = [value for key, value in path_different]

        for path in path_different.values():
            expected.add_path(path)

        self.assertTrue(
            nx.is_isomorphic(
                build_splicegraph(exons),
                expected
            )
        )


if __name__ == '__main__':
    unittest.main()
