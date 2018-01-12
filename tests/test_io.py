#!/usr/bin/env python3

"""
Tests for exfi.io
"""


from unittest import TestCase, main

import filecmp
import tempfile
import os

from exfi.io import \
    _coordinate_str_to_tuple, \
    _compute_segments, \
    _compute_links, \
    _compute_containments, \
    _compute_paths, \
    write_gfa1, \
    gfa1_to_exons, \
    gfa1_to_gapped_transcript

from tests.test_data import \
    EXONS_EMPTY_FN, EXONS_SIMPLE_FN, \
    EXONS_COMPLEX_FN, EXONS_COMPLEX_SOFT_FN, EXONS_COMPLEX_HARD_FN, \
    TRANSCRIPTOME_EMPTY, TRANSCRIPTOME_SIMPLE, TRANSCRIPTOME_COMPLEX, \
    GFA_EMPTY_FN, GFA_SIMPLE_FN, GFA_COMPLEX_FN, \
    GAPPED_EMPTY_FN, GAPPED_SIMPLE_FN, \
    GAPPED_COMPLEX_FN, GAPPED_COMPLEX_HARD_FN, GAPPED_COMPLEX_SOFT_FN, \
    SPLICE_GRAPH_EMPTY, SPLICE_GRAPH_SIMPLE, SPLICE_GRAPH_COMPLEX, \
    SEGMENTS_EMPTY, SEGMENTS_SIMPLE, SEGMENTS_COMPLEX, \
    LINKS_EMPTY, LINKS_SIMPLE, LINKS_COMPLEX, \
    CONTAINMENTS_EMPTY, CONTAINMENTS_SIMPLE, CONTAINMENTS_COMPLEX, \
    PATHS_EMPTY, PATHS_SIMPLE, PATHS_COMPLEX


class TestCoordinatesToVariables(TestCase):
    """_coordinate_str_to_tuple(coordinates):
    (string) -> (string, int, int)
    """
    def test_empty_string(self):
        """_coordinate_str_to_tuple: case empty"""
        with self.assertRaises(ValueError):
            _coordinate_str_to_tuple("")

    def test_incorrect_string(self):
        """_coordinate_str_to_tuple: case incorrect"""
        with self.assertRaises(IndexError):
            _coordinate_str_to_tuple("ENSDART00000161035:1")

    def test_correct_string(self):
        """_coordinate_str_to_tuple: case correct"""
        self.assertEqual(
            _coordinate_str_to_tuple("ENSDART00000161035:1-15"),
            ("ENSDART00000161035", 1, 15)
        )

    def test_messy_string(self):
        """_coordinate_str_to_tuple: messy case"""
        self.assertEqual(
            _coordinate_str_to_tuple("TRINITY_g14_c15_i5:1:2-15-12:1-15"),
            ("TRINITY_g14_c15_i5:1:2-15-12", 1, 15)
        )


class TestComputeSegments(TestCase):
    """Tests for _compute_segments"""

    def test_empty(self):
        """_compute_segments: empty case"""
        actual = list(_compute_segments(SPLICE_GRAPH_EMPTY, TRANSCRIPTOME_EMPTY))
        expected = SEGMENTS_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """_compute_segments: simple case"""
        actual = list(_compute_segments(SPLICE_GRAPH_SIMPLE, TRANSCRIPTOME_SIMPLE))
        expected = SEGMENTS_SIMPLE
        self.assertEqual(actual, expected)

    def test_complex(self):
        """_compute_segments: complex case"""
        actual = list(_compute_segments(SPLICE_GRAPH_COMPLEX, TRANSCRIPTOME_COMPLEX))
        expected = SEGMENTS_COMPLEX
        self.assertEqual(actual, expected)



class TestComputeLinks(TestCase):
    """Tests for _compute_links"""

    def test_empty(self):
        """_compute_links: empty case"""
        actual = list(_compute_links(SPLICE_GRAPH_EMPTY))
        expected = LINKS_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """_compute_links: simple case"""
        actual = list(_compute_links(SPLICE_GRAPH_SIMPLE))
        expected = LINKS_SIMPLE
        self.assertEqual(actual, expected)

    def test_coplex(self):
        """_compute_links: complex case"""
        actual = list(_compute_links(SPLICE_GRAPH_COMPLEX))
        expected = LINKS_COMPLEX
        self.assertEqual(actual, expected)



class TestComputeContainments(TestCase):
    """Tests for _compute_containments"""

    def test_empty(self):
        """_compute_containments: empty case"""
        actual = list(_compute_containments(SPLICE_GRAPH_EMPTY, TRANSCRIPTOME_EMPTY))
        expected = CONTAINMENTS_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """_compute_containments: simple case"""
        actual = list(_compute_containments(SPLICE_GRAPH_SIMPLE, TRANSCRIPTOME_SIMPLE))
        expected = CONTAINMENTS_SIMPLE
        self.assertEqual(actual, expected)

    def test_coplex(self):
        """_compute_containments: complex case"""
        actual = list(_compute_containments(SPLICE_GRAPH_COMPLEX, TRANSCRIPTOME_COMPLEX))
        expected = CONTAINMENTS_COMPLEX
        self.assertEqual(actual, expected)



class TestComputePaths(TestCase):
    """Tests for _compute_paths"""

    def test_empty(self):
        """_compute_paths: empty case"""
        actual = list(_compute_paths(SPLICE_GRAPH_EMPTY))
        expected = PATHS_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """_compute_paths: simple case"""
        actual = list(_compute_paths(SPLICE_GRAPH_SIMPLE))
        expected = PATHS_SIMPLE
        self.assertEqual(actual, expected)

    def test_complex(self):
        """_compute_paths: complex case"""
        actual = list(_compute_paths(SPLICE_GRAPH_COMPLEX))
        expected = PATHS_COMPLEX
        self.assertEqual(actual, expected)



class TestWriteGFA1(TestCase):
    """Tests for write_gfa1"""

    def test_empty(self):
        """write_gfa1: empty case"""
        tmp_file = tempfile.mkstemp()[1]
        write_gfa1(
            splice_graph=SPLICE_GRAPH_EMPTY,
            transcriptome_dict=TRANSCRIPTOME_EMPTY,
            filename=tmp_file
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            GFA_EMPTY_FN
        ))
        os.remove(tmp_file)

    def test_simple(self):
        """write_gfa1: simple case"""
        tmp_file = tempfile.mkstemp()[1]
        write_gfa1(
            splice_graph=SPLICE_GRAPH_SIMPLE,
            transcriptome_dict=TRANSCRIPTOME_SIMPLE,
            filename=tmp_file #tmp_file
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            GFA_SIMPLE_FN
        ))
        os.remove(tmp_file)

    def test_multiple(self):
        """write_gfa1: complex case"""
        tmp_file = tempfile.mkstemp()[1]
        write_gfa1(
            splice_graph=SPLICE_GRAPH_COMPLEX,
            transcriptome_dict=TRANSCRIPTOME_COMPLEX,
            filename=tmp_file
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            GFA_COMPLEX_FN
        ))
        os.remove(tmp_file)



class TestGFA1ToExons(TestCase):
    """Tests for gfa1_to_exons"""

    def test_empty(self):
        """gfa1_to_exons: empty case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=GFA_EMPTY_FN,
            fasta_out_fn=tmp_file,
            soft_mask_overlaps=False
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            EXONS_EMPTY_FN
        ))
        os.remove(tmp_file)


    def test_simple(self):
        """gfa1_to_exons: simple case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=GFA_SIMPLE_FN,
            fasta_out_fn=tmp_file,
            soft_mask_overlaps=False
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            EXONS_SIMPLE_FN
        ))
        os.remove(tmp_file)

    def test_multiple(self):
        """gfa1_to_exons: complex case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(gfa_in_fn=GFA_COMPLEX_FN, fasta_out_fn=tmp_file)
        self.assertTrue(filecmp.cmp(
            tmp_file, EXONS_COMPLEX_FN
        ))
        os.remove(tmp_file)

    def test_multiple_soft(self):
        """gfa1_to_exons: complex case and soft masking"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=GFA_COMPLEX_FN,
            fasta_out_fn=tmp_file,
            soft_mask_overlaps=True
        )
        self.assertTrue(
            filecmp.cmp(tmp_file, EXONS_COMPLEX_SOFT_FN)
        )
        os.remove(tmp_file)

    def test_multiple_hard(self):
        """gfa1_to_exons: complex case and hard masking"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=GFA_COMPLEX_FN,
            fasta_out_fn=tmp_file,
            hard_mask_overlaps=True
        )
        self.assertTrue(
            filecmp.cmp(tmp_file, EXONS_COMPLEX_HARD_FN)
        )
        os.remove(tmp_file)



class TestGFA1ToGappedTranscript(TestCase):
    """Tests for gfa1_to_gapped_transcript"""

    def test_empty(self):
        """gfa1_to_gapped_transcript: empty case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(gfa_in=GFA_EMPTY_FN, fasta_out=tmp_file)
        self.assertTrue(
            filecmp.cmp(tmp_file, GAPPED_EMPTY_FN)
        )
        os.remove(tmp_file)

    def test_simple(self):
        """gfa1_to_gapped_transcript: simple case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(gfa_in=GFA_SIMPLE_FN, fasta_out=tmp_file)
        self.assertTrue(
            filecmp.cmp(tmp_file, GAPPED_SIMPLE_FN)
        )
        os.remove(tmp_file)

    def test_multiple(self):
        """gfa1_to_gapped_transcript: complex case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(gfa_in=GFA_COMPLEX_FN, fasta_out=tmp_file)
        self.assertTrue(filecmp.cmp(
            tmp_file,
            GAPPED_COMPLEX_FN
        ))
        os.remove(tmp_file)

    def test_multiple_soft(self):
        """gfa1_to_gapped_transcript: complex case and soft masking"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(gfa_in=GFA_COMPLEX_FN, fasta_out=tmp_file,
                                  soft_mask_overlaps=True)
        self.assertTrue(
            filecmp.cmp(tmp_file, GAPPED_COMPLEX_SOFT_FN)
        )
        os.remove(tmp_file)


    def test_multiple_hard(self):
        """gfa1_to_gapped_transcript: complex case and hard masking"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(gfa_in=GFA_COMPLEX_FN, fasta_out=tmp_file,
                                  hard_mask_overlaps=True)
        self.assertTrue(
            filecmp.cmp(tmp_file, GAPPED_COMPLEX_HARD_FN)
        )
        os.remove(tmp_file)

if __name__ == '__main__':
    main()
