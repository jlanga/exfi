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
    exons_empty, exons_simple, exons_complex, \
    exons_complex_soft, exons_complex_hard, \
    transcriptome_empty, transcriptome_simple, transcriptome_complex, \
    gfa_empty, gfa_simple, gfa_complex, \
    gapped_empty, gapped_simple, gapped_complex, \
    gapped_complex_hard, gapped_complex_soft, \
    splice_graph_empty, splice_graph_simple, splice_graph_complex, \
    segments_empty, segments_simple, segments_complex, \
    links_empty, links_simple, links_complex, \
    containments_empty, containments_simple, containments_complex, \
    paths_empty, paths_simple, paths_complex



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
        actual = list(_compute_segments(splice_graph_empty, transcriptome_empty))
        expected = segments_empty
        self.assertEqual(actual, expected)

    def test_simple(self):
        """_compute_segments: simple case"""
        actual = list(_compute_segments(splice_graph_simple, transcriptome_simple))
        expected = segments_simple
        self.assertEqual(actual, expected)

    def test_complex(self):
        """_compute_segments: complex case"""
        actual = list(_compute_segments(splice_graph_complex, transcriptome_complex))
        expected = segments_complex
        self.assertEqual(actual, expected)



class TestComputeLinks(TestCase):
    """Tests for _compute_links"""

    def test_empty(self):
        """_compute_links: empty case"""
        actual = list(_compute_links(splice_graph_empty))
        expected = links_empty
        self.assertEqual(actual, expected)

    def test_simple(self):
        """_compute_links: simple case"""
        actual = list(_compute_links(splice_graph_simple))
        expected = links_simple
        self.assertEqual(actual, expected)

    def test_coplex(self):
        """_compute_links: complex case"""
        actual = list(_compute_links(splice_graph_complex))
        expected = links_complex
        self.assertEqual(actual, expected)



class TestComputeContainments(TestCase):
    """Tests for _compute_containments"""

    def test_empty(self):
        """_compute_containments: empty case"""
        actual = list(_compute_containments(splice_graph_empty, transcriptome_empty))
        expected = containments_empty
        self.assertEqual(actual, expected)

    def test_simple(self):
        """_compute_containments: simple case"""
        actual = list(_compute_containments(splice_graph_simple, transcriptome_simple))
        expected = containments_simple
        self.assertEqual(actual, expected)

    def test_coplex(self):
        """_compute_containments: complex case"""
        actual = list(_compute_containments(splice_graph_complex, transcriptome_complex))
        expected = containments_complex
        self.assertEqual(actual, expected)



class TestComputePaths(TestCase):
    """Tests for _compute_paths"""

    def test_empty(self):
        """_compute_paths: empty case"""
        actual = list(_compute_paths(splice_graph_empty))
        expected = paths_empty
        self.assertEqual(actual, expected)

    def test_simple(self):
        """_compute_paths: simple case"""
        actual = list(_compute_paths(splice_graph_simple))
        expected = paths_simple
        self.assertEqual(actual, expected)

    def test_complex(self):
        """_compute_paths: complex case"""
        actual = list(_compute_paths(splice_graph_complex))
        expected = paths_complex
        self.assertEqual(actual, expected)



class TestWriteGFA1(TestCase):
    """Tests for write_gfa1"""

    def test_empty(self):
        """write_gfa1: empty case"""
        tmp_file = tempfile.mkstemp()[1]
        write_gfa1(
            splice_graph=splice_graph_empty,
            transcriptome_dict=transcriptome_empty,
            filename=tmp_file
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            gfa_empty
        ))
        os.remove(tmp_file)

    def test_simple(self):
        """write_gfa1: simple case"""
        tmp_file = tempfile.mkstemp()[1]
        write_gfa1(
            splice_graph=splice_graph_simple,
            transcriptome_dict=transcriptome_simple,
            filename=tmp_file #tmp_file
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            gfa_simple
        ))
        os.remove(tmp_file)

    def test_multiple(self):
        """write_gfa1: complex case"""
        tmp_file = tempfile.mkstemp()[1]
        write_gfa1(
            splice_graph=splice_graph_complex,
            transcriptome_dict=transcriptome_complex,
            filename=tmp_file
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            gfa_complex
        ))
        os.remove(tmp_file)



class TestGFA1ToExons(TestCase):
    """Tests for gfa1_to_exons"""

    def test_empty(self):
        """gfa1_to_exons: empty case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=gfa_empty,
            fasta_out_fn=tmp_file,
            soft_mask_overlaps=False
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            exons_empty
        ))
        os.remove(tmp_file)


    def test_simple(self):
        """gfa1_to_exons: simple case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=gfa_simple,
            fasta_out_fn=tmp_file,
            soft_mask_overlaps=False
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            exons_simple
        ))
        os.remove(tmp_file)

    def test_multiple(self):
        """gfa1_to_exons: complex case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(gfa_in_fn=gfa_complex, fasta_out_fn=tmp_file)
        self.assertTrue(
            filecmp.cmp(tmp_file, exons_complex)
        )
        os.remove(tmp_file)

    def test_multiple_soft(self):
        """gfa1_to_exons: complex case and soft masking"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=gfa_complex,
            fasta_out_fn=tmp_file,
            soft_mask_overlaps=True
        )
        self.assertTrue(
            filecmp.cmp(tmp_file, exons_complex_soft)
        )
        os.remove(tmp_file)

    def test_multiple_hard(self):
        """gfa1_to_exons: complex case and hard masking"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=gfa_complex,
            fasta_out_fn=tmp_file,
            hard_mask_overlaps=True
        )
        self.assertTrue(
            filecmp.cmp(tmp_file, exons_complex_hard)
        )
        os.remove(tmp_file)



class TestGFA1ToGappedTranscript(TestCase):
    """Tests for gfa1_to_gapped_transcript"""

    def test_empty(self):
        """gfa1_to_gapped_transcript: empty case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(gfa_in=gfa_empty, fasta_out=tmp_file)
        self.assertTrue(
            filecmp.cmp(tmp_file, gapped_empty)
        )
        os.remove(tmp_file)

    def test_simple(self):
        """gfa1_to_gapped_transcript: simple case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(gfa_in=gfa_simple, fasta_out=tmp_file)
        self.assertTrue(
            filecmp.cmp(tmp_file, gapped_simple)
        )
        os.remove(tmp_file)

    def test_multiple(self):
        """gfa1_to_gapped_transcript: complex case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(gfa_in=gfa_complex, fasta_out=tmp_file)
        self.assertTrue(filecmp.cmp(
            tmp_file,
            gapped_complex
        ))
        os.remove(tmp_file)

    def test_multiple_soft(self):
        """gfa1_to_gapped_transcript: complex case and soft masking"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(gfa_in=gfa_complex, fasta_out=tmp_file,
                                  soft_mask_overlaps=True)
        self.assertTrue(
            filecmp.cmp(tmp_file, gapped_complex_soft)
        )
        os.remove(tmp_file)


    def test_multiple_hard(self):
        """gfa1_to_gapped_transcript: complex case and hard masking"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(gfa_in=gfa_complex, fasta_out=tmp_file,
                                  hard_mask_overlaps=True)
        self.assertTrue(
            filecmp.cmp(tmp_file, gapped_complex_hard)
        )
        os.remove(tmp_file)

if __name__ == '__main__':
    main()
