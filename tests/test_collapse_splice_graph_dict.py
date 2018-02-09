#!/usr/bin/env python3

"""
Tests for exfi.collapse_splice_graph
"""

from unittest import TestCase, main

import networkx as nx

from exfi.collapse_splice_graph_dict import \
    _compute_seq2node, \
    _compute_old2new, \
    _compute_new_node2coord, \
    _compute_new_link2overlap, \
    collapse_splice_graph_dict

from tests.custom_assertions import \
    CustomAssertions

from tests.test_data import \
    NODE2COORDS_EMPTY, NODE2COORDS_SIMPLE, NODE2COORDS_COMPLEX, \
    OVERLAPS_EMPTY, OVERLAPS_SIMPLE, OVERLAPS_COMPLEX, \
    TRANSCRIPTOME_EMPTY_DICT, TRANSCRIPTOME_SIMPLE_DICT, TRANSCRIPTOME_COMPLEX_DICT, \
    SPLICE_GRAPH_EMPTY_DICT, SPLICE_GRAPH_SIMPLE_DICT, SPLICE_GRAPH_COMPLEX_DICT


SEQ2NODE_EMPTY = {}
SEQ2NODE_SIMPLE = {
    'TGCACGGGTTTATTGTTCACAAAGAGATCGACAATGTGCGCAACTAAAATAAACATAGTACATTTTGATTATACACGAACTTAAACTAAAGTCC'
    'AATCACACCTCCGCCCCGTTTCCACAGCAGCCTGTCAGGGTGGAGGAAAAGCGCGGCGGTCATGTGAGGCTCGAGCATCTCTCTCTCTCTCTCT'
    'CTCTCTCTCTCTACAGAATGATAGAGGGAGCTCGTGAATCACATCATAGTCGTCCTCCCCTCATTCGTCCTCTCCAGCAGACACCGAAAAACTG'
    'CGTTCATGCCAAAATGGGATGTGGAAATTCCTCCGCCACGAGCA':
        ('ENSDART00000161035.1:0-326',)
}
SEQ2NODE_COMPLEX = {
    "TGCACGGGTTTATTGTTCACAAAGAGATCGACAATGTGCGCAACTAAAATAAACATAGTACATTTTGATTATACACGAACTTAAACTAAAGTCC"
    "AATCACACCTCCGCCCCGTTTCCACAGCAGCCTGTCAGGGTGGAGGAAAAGCGCGGCGGTCATGTGAGGCTCGAGCATCTCTCTCTCTCTCTCT"
    "CTCTCTCTCTCTACAGAATGATAGAGGGAGCTCGTGAATCACATCATAGTCGTCCTCCCCTCATTCGTCCTCTCCAGCAGACACCGAAAAACTG"
    "CGTTCATGCCAAAATGGGATGTGGAAATTCCTCCGCCACGAGCA":
        ("ENSDART00000161035.1:0-326", ),
    "AGGAACTACGGTGGAGTGTATGTGGGTCTTCCTGCTGATCTGACTGCAGTCGCTGCCAGTCAGTCCAAATCAACA":
        ("ENSDART00000161035.1:397-472", ),
    "AGTCAACAGATGTTTATTGCAGACCTTCAGATAAAACAACATAGAA":
        ("ENSDART00000161035.1:477-523", ),
    "TGGAGCTGAAGCCGAGTATCTTGGTATTGGACTGGAACAGAAATCCAGCAAAAACTTTAAGGGAAATCACTTTCATTTCATGATCGAAAAACTC"
    "CCGCAGATCATAAAAGAGTGGAAGGAAG":
        ("ENSDART00000165342.1:5-127", ),
    "AGGACCTGTAGTAGAAACAAAACTAGGATCTCTGAGAGGTGCCTTCTTGACTGTGAAGGGCAAGGACACAATAGTCAATAGTTATCTAGGTGTG"
    "CCGTTCGCCAAGCCGCCTGTAGGACCCCTGAGACTTGCTCGACCACAGGCTGCAGAGAAATGGCAAGGAGTTAGAGATGCCACCA":
        ("ENSDART00000165342.1:125-304", ),
    "GTGCCTCCAGGAAAGGCAAATGACTGTAACTGAACTGGAGTTTCTATCGATGGATGTGGAGGTTCCTGAGGTCTCGGAGGATTGCCTGTATCTT"
    "AACATCTACACCCCAGTTAAACCTGGACAAGGAGACAAGAAGTTACCAG":
        ("ENSDART00000165342.1:317-460", ),
    "GTCATGGTTTGGATTCATGGTGGAGGACTCTCTCTTGGATCGGCTTCAATGTATGATGGCTCTGTTCTGGCTGCGTATCAGGATGTGGTCGTGG"
    "TGCTCATTCAGTACAGATTGGGTCTTCTGGGGTTCTTAA":
        ("ENSDART00000165342.1:459-592", ),
    "AGCACCGGAGACGAGCATGCGCCAGGAAACTATGGTTTTCTGGATCAAGTAGCTGCCCT":
        ("ENSDART00000165342.1:591-650", ),
    "GCCCTTCAGTGGGTTCAGGAGAACATCCACAGCTTCGGTGGAGATCCTGGATCAGTGACCATCTTTGGAGAGTCTGCTGGAGGAATCAGTGTAT"
    "CCACGCT":
        ("ENSDART00000165342.1:645-746", ),
    "GATTCTTTCCCCGCTGGCGTCTGGACTGTTTCATCGCGCCATTGCAGAAAGTGGAACTGCCTTCTGGGATGGTTTAGTCATGGCTGATCCTTTT"
    "CAGAGAGCCCA":
        ("ENSDART00000165342.1:746-851", ),
    "TGCAGCCAAACAATGCAACTGTGACAGCAGCA":
        ("ENSDART00000165342.1:854-886", ),
    "TGTCGACTGCATTATGCACTGGTCTGAAGAGGAGGCTCTGGAATGTGCTAAAAA":
        ("ENSDART00000165342.1:899-953", ),
    "CGTTGCTGTAGATTCTTATTTCCTTCCCAAACCCATCGAGGAGATTGTTGAGAAACAAGAGTTTAGTAAAGTTCCTCTCATCAACGGCATTAAC"
    "AATGATGAGTTTGGCTTCTTGTTGGCTGA":
        ("ENSDART00000165342.1:974-1097", ),
    "TATTTCTTGGGTCCTGAATGGATGAATGGGTTGAAAAGAGAGCAAATCGCTGAAGCCTTGACGCTCACATATCCTGA":
        ("ENSDART00000165342.1:1098-1175", ),
    "CCCAAGGATCGATGGATCATTGATCTGGTGGCGAAGGAATATCTGGGCGACACACACGACCCCATTGAAATCCGTGAAGTTTATCGGGAGATGA"
    "TGGGAGACGTGCTGTTTAACATCCCTGCCCTGCAACTGGCAAAACACCACAGCG":
        ("ENSDART00000165342.1:1176-1324", )
}

OLD2NEW_EMPTY = {}
OLD2NEW_SIMPLE = {
    'ENSDART00000161035.1:0-326': 'exon_00000000'
}
OLD2NEW_COMPLEX = {
    'ENSDART00000161035.1:0-326': 'exon_00000000',
    'ENSDART00000161035.1:397-472': 'exon_00000001',
    'ENSDART00000161035.1:477-523': 'exon_00000002',
    'ENSDART00000165342.1:5-127': 'exon_00000003',
    'ENSDART00000165342.1:125-304': 'exon_00000004',
    'ENSDART00000165342.1:317-460': 'exon_00000005',
    'ENSDART00000165342.1:459-592': 'exon_00000006',
    'ENSDART00000165342.1:591-650': 'exon_00000007',
    'ENSDART00000165342.1:645-746': 'exon_00000008',
    'ENSDART00000165342.1:746-851': 'exon_00000009',
    'ENSDART00000165342.1:854-886': 'exon_00000010',
    'ENSDART00000165342.1:899-953': 'exon_00000011',
    'ENSDART00000165342.1:974-1097': 'exon_00000012',
    'ENSDART00000165342.1:1098-1175': 'exon_00000013',
    'ENSDART00000165342.1:1176-1324': 'exon_00000014'
}

NEW_NODE2COORD_EMPTY = {}
NEW_NODE2COORD_SIMPLE = {
    'exon_00000000': (('ENSDART00000161035.1', 0, 326),)
}
NEW_NODE2COORD_COMPLEX = {
    'exon_00000000': (('ENSDART00000161035.1', 0, 326),),
    'exon_00000001': (('ENSDART00000161035.1', 397, 472),),
    'exon_00000002': (('ENSDART00000161035.1', 477, 523),),
    'exon_00000003': (('ENSDART00000165342.1', 5, 127),),
    'exon_00000004': (('ENSDART00000165342.1', 125, 304),),
    'exon_00000005': (('ENSDART00000165342.1', 317, 460),),
    'exon_00000006': (('ENSDART00000165342.1', 459, 592),),
    'exon_00000007': (('ENSDART00000165342.1', 591, 650),),
    'exon_00000008': (('ENSDART00000165342.1', 645, 746),),
    'exon_00000009': (('ENSDART00000165342.1', 746, 851),),
    'exon_00000010': (('ENSDART00000165342.1', 854, 886),),
    'exon_00000011': (('ENSDART00000165342.1', 899, 953),),
    'exon_00000012': (('ENSDART00000165342.1', 974, 1097),),
    'exon_00000013': (('ENSDART00000165342.1', 1098, 1175),),
    'exon_00000014': (('ENSDART00000165342.1', 1176, 1324),)
}

LINK2OVERLAP_EMPTY = OVERLAPS_EMPTY
LINK2OVERLAP_SIMPLE = OVERLAPS_SIMPLE
LINK2OVERLAP_COMPLEX = OVERLAPS_COMPLEX

NEW_LINK2OVERLAP_EMPTY = {}
NEW_LINK2OVERLAP_SIMPLE = {}
NEW_LINK2OVERLAP_COMPLEX = {
    ('exon_00000000', 'exon_00000001'): -71,
    ('exon_00000001', 'exon_00000002'): -5,
    ('exon_00000003', 'exon_00000004'): 2,
    ('exon_00000004', 'exon_00000005'): -13,
    ('exon_00000005', 'exon_00000006'): 1,
    ('exon_00000006', 'exon_00000007'): 1,
    ('exon_00000007', 'exon_00000008'): 5,
    ('exon_00000008', 'exon_00000009'): 0,
    ('exon_00000009', 'exon_00000010'): -3,
    ('exon_00000010', 'exon_00000011'): -13,
    ('exon_00000011', 'exon_00000012'): -21,
    ('exon_00000012', 'exon_00000013'): -1,
    ('exon_00000013', 'exon_00000014'): -1
}


COLLAPSED_EMPTY = nx.DiGraph()

COLLAPSED_SIMPLE = nx.DiGraph()
COLLAPSED_SIMPLE.add_nodes_from(NEW_NODE2COORD_SIMPLE.keys())
COLLAPSED_SIMPLE.add_edges_from(NEW_LINK2OVERLAP_SIMPLE.keys())
nx.set_node_attributes(G=COLLAPSED_SIMPLE, name="coordinates", values=NEW_NODE2COORD_SIMPLE)
nx.set_edge_attributes(G=COLLAPSED_SIMPLE, name="overlaps", values=NEW_LINK2OVERLAP_SIMPLE)

COLLAPSED_COMPLEX = nx.DiGraph()
COLLAPSED_COMPLEX.add_nodes_from(NEW_NODE2COORD_COMPLEX.keys())
COLLAPSED_COMPLEX.add_edges_from(NEW_LINK2OVERLAP_COMPLEX.keys())
nx.set_node_attributes(G=COLLAPSED_COMPLEX, name="coordinates", values=NEW_NODE2COORD_COMPLEX)
nx.set_edge_attributes(G=COLLAPSED_COMPLEX, name="overlaps", values=NEW_LINK2OVERLAP_COMPLEX)



class TestComputeSeq2Node(TestCase):
    """Tests for exfi.collapse_splice_graph._compute_seq2node"""

    def test_empty(self):
        """exfi.collapse_splice_graph._compute_seq2node: empty case"""
        actual = _compute_seq2node(NODE2COORDS_EMPTY, TRANSCRIPTOME_EMPTY_DICT)
        expected = SEQ2NODE_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """exfi.collapse_splice_graph._compute_seq2node: simple case"""
        actual = _compute_seq2node(NODE2COORDS_SIMPLE, TRANSCRIPTOME_SIMPLE_DICT)
        expected = SEQ2NODE_SIMPLE
        self.assertEqual(actual, expected)

    def test_complex(self):
        """exfi.collapse_splice_graph._compute_seq2node: complex case"""
        actual = _compute_seq2node(NODE2COORDS_COMPLEX, TRANSCRIPTOME_COMPLEX_DICT)
        expected = SEQ2NODE_COMPLEX
        self.assertEqual(actual, expected)



class TestComputeOld2New(TestCase):
    """Tests for exfi.collapse_splice_graph._compute_old2new"""

    def test_empty(self):
        """exfi.collapse_splice_graph._compute_old2new: empty case"""
        actual = _compute_old2new(SEQ2NODE_EMPTY)
        expected = OLD2NEW_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """exfi.collapse_splice_graph._compute_seq2node: simple case"""
        actual = _compute_old2new(SEQ2NODE_SIMPLE)
        expected = OLD2NEW_SIMPLE
        self.assertEqual(actual, expected)

    def test_complex(self):
        """exfi.collapse_splice_graph._compute_seq2node: complex case"""
        actual = _compute_old2new(SEQ2NODE_COMPLEX)
        expected = OLD2NEW_COMPLEX
        self.assertEqual(actual, expected)



class TestComputeNewNode2Coord(TestCase):
    """Tests for exfi.collapse_splice_graph._compute_new_node2coord"""

    def test_empty(self):
        """exfi.collapse_splice_graph._compute_new_node2coord: empty case"""
        actual = _compute_new_node2coord(OLD2NEW_EMPTY, NODE2COORDS_EMPTY)
        expected = NEW_NODE2COORD_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """exfi.collapse_splice_graph._compute_new_node2coord: simple case"""
        actual = _compute_new_node2coord(OLD2NEW_SIMPLE, NODE2COORDS_SIMPLE)
        expected = NEW_NODE2COORD_SIMPLE
        self.assertEqual(actual, expected)

    def test_complex(self):
        """exfi.collapse_splice_graph._compute_new_node2coord: complex case"""
        actual = _compute_new_node2coord(OLD2NEW_COMPLEX, NODE2COORDS_COMPLEX)
        expected = NEW_NODE2COORD_COMPLEX
        self.assertEqual(actual, expected)



class TestComputeNewLink2Overlap(TestCase):
    """Tests for exfi.collapse_splice_graph._compute_new_link2overlap"""

    def test_empty(self):
        """exfi.collapse_splice_graph._compute_new_link2overlap: empty case"""
        actual = _compute_new_link2overlap(OLD2NEW_EMPTY, LINK2OVERLAP_EMPTY)
        expected = NEW_LINK2OVERLAP_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """exfi.collapse_splice_graph._compute_new_link2overlap: simple case"""
        actual = _compute_new_link2overlap(OLD2NEW_SIMPLE, LINK2OVERLAP_SIMPLE)
        expected = NEW_LINK2OVERLAP_SIMPLE
        self.assertEqual(actual, expected)

    def test_complex(self):
        """exfi.collapse_splice_graph._compute_new_link2overlap: complex case"""
        actual = _compute_new_link2overlap(OLD2NEW_COMPLEX, LINK2OVERLAP_COMPLEX)
        expected = NEW_LINK2OVERLAP_COMPLEX
        self.assertEqual(actual, expected)




class TestCollapseSpliceGraph(TestCase, CustomAssertions):
    """Tests for exfi.collapse_splice_graph_dict.collapse_splice_graph_dict"""

    def test_empty(self):
        """exfi.collapse_splice_graph.collapse_splice_graph: empty case"""
        actual = collapse_splice_graph_dict(SPLICE_GRAPH_EMPTY_DICT, TRANSCRIPTOME_EMPTY_DICT)
        expected = COLLAPSED_EMPTY
        self.assertEqualSpliceGraphs(actual, expected)

    def test_simple(self):
        """exfi.collapse_splice_graph.collapse_splice_graph: simple case"""
        actual = collapse_splice_graph_dict(SPLICE_GRAPH_SIMPLE_DICT, TRANSCRIPTOME_SIMPLE_DICT)
        expected = COLLAPSED_SIMPLE
        self.assertEqualSpliceGraphs(actual, expected)

    def test_complex(self):
        """exfi.collapse_splice_graph.collapse_splice_graph: complex case"""
        actual = collapse_splice_graph_dict(SPLICE_GRAPH_COMPLEX_DICT, TRANSCRIPTOME_COMPLEX_DICT)
        expected = COLLAPSED_COMPLEX
        self.assertEqualSpliceGraphs(actual, expected)



if __name__ == "__main__":
    main()
