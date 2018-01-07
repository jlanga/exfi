#!/usr/bin/env python3

import unittest

from exfi.io import \
    _compute_segments, \
    _compute_links, \
    _compute_containments, \
    _compute_paths, \
    write_gfa1, \
    gfa1_to_exons, \
    gfa1_to_gapped_transcript

import filecmp
import tempfile
import os

from tests.test_data import *


class TestComputeSegments(unittest.TestCase):

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



class TestComputeLinks(unittest.TestCase):
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



class TestComputeContainments(unittest.TestCase):
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



class TestComputePaths(unittest.TestCase):
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



class TestWriteGFA1(unittest.TestCase):

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
            empty_gfa
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
            simple_gfa
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
            complex_gfa
        ))
        os.remove(tmp_file)



class TestGFA1ToExons(unittest.TestCase):

    def test_empty(self):
        """gfa1_to_exons: empty case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=empty_gfa,
            fasta_out_fn=tmp_file,
            soft_mask_overlaps=False
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            empty_exons
        ))
        os.remove(tmp_file)


    def test_simple(self):
        """gfa1_to_exons: simple case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=simple_gfa,
            fasta_out_fn=tmp_file,
            soft_mask_overlaps=False
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            simple_exons
        ))
        os.remove(tmp_file)

    def test_multiple(self):
        """gfa1_to_exons: complex case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(gfa_in_fn=complex_gfa, fasta_out_fn=tmp_file)
        self.assertTrue(
            filecmp.cmp(tmp_file, complex_exons)
        )
        os.remove(tmp_file)

    def test_multiple_soft(self):
        """gfa1_to_exons: complex case and soft masking"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=complex_gfa,
            fasta_out_fn=tmp_file,
            soft_mask_overlaps=True
        )
        self.assertTrue(
            filecmp.cmp(tmp_file, complex_exons_soft)
        )
        os.remove(tmp_file)

    def test_multiple_hard(self):
        """gfa1_to_exons: complex case and hard masking"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=complex_gfa,
            fasta_out_fn=tmp_file,
            hard_mask_overlaps=True
        )
        self.assertTrue(
            filecmp.cmp(tmp_file, complex_exons_hard)
        )
        os.remove(tmp_file)



class TestGFA1ToGappedTranscript(unittest.TestCase):

    def test_empty(self):
        """gfa1_to_gapped_transcript: empty case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(gfa_in=empty_gfa, fasta_out=tmp_file)
        self.assertTrue(
            filecmp.cmp(tmp_file, empty_gapped)
        )
        os.remove(tmp_file)

    def test_simple(self):
        """gfa1_to_gapped_transcript: simple case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(gfa_in=simple_gfa, fasta_out=tmp_file)
        self.assertTrue(
            filecmp.cmp(tmp_file, simple_gapped)
        )
        os.remove(tmp_file)

    def test_multiple(self):
        """gfa1_to_gapped_transcript: complex case"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(gfa_in=complex_gfa, fasta_out=tmp_file)
        self.assertTrue(filecmp.cmp(
            tmp_file,
            complex_gapped
        ))
        os.remove(tmp_file)

    def test_multiple_soft(self):
        """gfa1_to_gapped_transcript: complex case and soft masking"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(gfa_in=complex_gfa, fasta_out=tmp_file,
            soft_mask_overlaps=True
        )
        self.assertTrue(
            filecmp.cmp(tmp_file, complex_gapped_soft)
        )
        os.remove(tmp_file)


    def test_multiple_hard(self):
        """gfa1_to_gapped_transcript: complex case and hard masking"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(gfa_in=complex_gfa, fasta_out=tmp_file,
            hard_mask_overlaps=True
        )
        self.assertTrue(
            filecmp.cmp(tmp_file, complex_gapped_hard)
        )
        os.remove(tmp_file)

if __name__ == '__main__':
    unittest.main()
