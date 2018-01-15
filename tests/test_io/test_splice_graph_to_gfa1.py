#!/usr/bin/env python3

"""
Tests for exfi.io.gfa1_to_exons
"""

from unittest import TestCase, main

import filecmp
import tempfile
import os

from exfi.io.splice_graph_to_gfa1 import \
    _compute_segments, \
    _compute_links, \
    _compute_containments, \
    _compute_paths, \
    splice_graph_to_gfa1

from tests.test_data import \
    TRANSCRIPTOME_EMPTY, TRANSCRIPTOME_SIMPLE, TRANSCRIPTOME_COMPLEX, \
    GFA_EMPTY_FN, GFA_SIMPLE_FN, GFA_COMPLEX_FN, \
    SPLICE_GRAPH_EMPTY, SPLICE_GRAPH_SIMPLE, SPLICE_GRAPH_COMPLEX, \
    SEGMENTS_EMPTY, SEGMENTS_SIMPLE, SEGMENTS_COMPLEX, \
    LINKS_EMPTY, LINKS_SIMPLE, LINKS_COMPLEX, \
    CONTAINMENTS_EMPTY, CONTAINMENTS_SIMPLE, CONTAINMENTS_COMPLEX, \
    PATHS_EMPTY, PATHS_SIMPLE, PATHS_COMPLEX



class TestComputeSegments(TestCase):
    """Tests for _compute_segments"""

    def test_empty(self):
        """exfi.io.splice_graph_to_gfa1._compute_segments: empty case"""
        actual = list(_compute_segments(SPLICE_GRAPH_EMPTY, TRANSCRIPTOME_EMPTY))
        expected = SEGMENTS_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """exfi.io.splice_graph_to_gfa1._compute_segments: simple case"""
        actual = list(_compute_segments(SPLICE_GRAPH_SIMPLE, TRANSCRIPTOME_SIMPLE))
        expected = SEGMENTS_SIMPLE
        self.assertEqual(actual, expected)

    def test_complex(self):
        """exfi.io.splice_graph_to_gfa1._compute_segments: complex case"""
        actual = list(_compute_segments(SPLICE_GRAPH_COMPLEX, TRANSCRIPTOME_COMPLEX))
        expected = SEGMENTS_COMPLEX
        self.assertEqual(actual, expected)



class TestComputeLinks(TestCase):
    """Tests for _compute_links"""

    def test_empty(self):
        """exfi.io.splice_graph_to_gfa1._compute_links: empty case"""
        actual = list(_compute_links(SPLICE_GRAPH_EMPTY))
        expected = LINKS_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """exfi.io.splice_graph_to_gfa1._compute_links: simple case"""
        actual = list(_compute_links(SPLICE_GRAPH_SIMPLE))
        expected = LINKS_SIMPLE
        self.assertEqual(actual, expected)

    def test_coplex(self):
        """exfi.io.splice_graph_to_gfa1._compute_links: complex case"""
        actual = list(_compute_links(SPLICE_GRAPH_COMPLEX))
        expected = LINKS_COMPLEX
        self.assertEqual(actual, expected)



class TestComputeContainments(TestCase):
    """Tests for _compute_containments"""

    def test_empty(self):
        """exfi.io.splice_graph_to_gfa1._compute_containments: empty case"""
        actual = list(_compute_containments(SPLICE_GRAPH_EMPTY, TRANSCRIPTOME_EMPTY))
        expected = CONTAINMENTS_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """exfi.io.splice_graph_to_gfa1._compute_containments: simple case"""
        actual = list(_compute_containments(SPLICE_GRAPH_SIMPLE, TRANSCRIPTOME_SIMPLE))
        expected = CONTAINMENTS_SIMPLE
        self.assertEqual(actual, expected)

    def test_coplex(self):
        """exfi.io.splice_graph_to_gfa1._compute_containments: complex case"""
        actual = list(_compute_containments(SPLICE_GRAPH_COMPLEX, TRANSCRIPTOME_COMPLEX))
        expected = CONTAINMENTS_COMPLEX
        self.assertEqual(actual, expected)



class TestComputePaths(TestCase):
    """Tests for exfi.io.splice_graph_to_gfa1._compute_paths"""

    def test_empty(self):
        """exfi.io.splice_graph_to_gfa1._compute_paths: empty case"""
        actual = list(_compute_paths(SPLICE_GRAPH_EMPTY))
        expected = PATHS_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """exfi.io.splice_graph_to_gfa1._compute_paths: simple case"""
        actual = list(_compute_paths(SPLICE_GRAPH_SIMPLE))
        expected = PATHS_SIMPLE
        self.assertEqual(actual, expected)

    def test_complex(self):
        """exfi.io.splice_graph_to_gfa1._compute_paths: complex case"""
        actual = list(_compute_paths(SPLICE_GRAPH_COMPLEX))
        expected = PATHS_COMPLEX
        self.assertEqual(actual, expected)


class TestSpliceGraphToGFA1(TestCase):
    """Tests for exfi.io.splice_graph_to_gfa1.splice_graph_to_gfa1"""

    def test_empty(self):
        """exfi.io.splice_graph_to_gfa1.splice_graph_to_gfa1: empty case"""
        tmp_file = tempfile.mkstemp()[1]
        splice_graph_to_gfa1(
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
        """exfi.io.splice_graph_to_gfa1.splice_graph_to_gfa1: simple case"""
        tmp_file = tempfile.mkstemp()[1]
        splice_graph_to_gfa1(
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
        """exfi.io.splice_graph_to_gfa1.splice_graph_to_gfa1: complex case"""
        tmp_file = tempfile.mkstemp()[1]
        splice_graph_to_gfa1(
            splice_graph=SPLICE_GRAPH_COMPLEX,
            transcriptome_dict=TRANSCRIPTOME_COMPLEX,
            filename=tmp_file
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            GFA_COMPLEX_FN
        ))
        os.remove(tmp_file)

if __name__ == "__main__":
    main()