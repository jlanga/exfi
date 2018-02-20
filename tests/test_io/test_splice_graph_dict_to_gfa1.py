#!/usr/bin/env python3

"""
Tests for exfi.io.gfa1_to_exons
"""

from unittest import TestCase, main

import filecmp
import tempfile
import os

from exfi.io.splice_graph_dict_to_gfa1 import \
    _compute_segments, \
    _compute_links, \
    _compute_containments, \
    _compute_paths, \
    splice_graph_dict_to_gfa1

from tests.test_data import \
    TRANSCRIPTOME_EMPTY_DICT, TRANSCRIPTOME_SIMPLE_DICT, TRANSCRIPTOME_COMPLEX_DICT, \
    GFA_EMPTY_FN, GFA_SIMPLE_FN, GFA_COMPLEX_FN, \
    SPLICE_GRAPH_EMPTY_DICT, SPLICE_GRAPH_SIMPLE_DICT, SPLICE_GRAPH_COMPLEX_DICT, \
    SEGMENTS_EMPTY, SEGMENTS_SIMPLE, SEGMENTS_COMPLEX, \
    LINKS_EMPTY, LINKS_SIMPLE, LINKS_COMPLEX, \
    CONTAINMENTS_EMPTY, CONTAINMENTS_SIMPLE, CONTAINMENTS_COMPLEX, \
    PATHS_EMPTY, PATHS_SIMPLE, PATHS_COMPLEX



class TestComputeSegments(TestCase):
    """Tests for _compute_segments"""

    def test_empty(self):
        """exfi.io.splice_graph_dict_to_gfa1._compute_segments: empty case"""
        actual = list(_compute_segments(SPLICE_GRAPH_EMPTY_DICT, TRANSCRIPTOME_EMPTY_DICT))
        expected = SEGMENTS_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """exfi.io.splice_graph_dict_to_gfa1._compute_segments: simple case"""
        actual = list(_compute_segments(SPLICE_GRAPH_SIMPLE_DICT, TRANSCRIPTOME_SIMPLE_DICT))
        expected = SEGMENTS_SIMPLE
        self.assertEqual(actual, expected)

    def test_complex(self):
        """exfi.io.splice_graph_dict_to_gfa1._compute_segments: complex case"""
        actual = list(_compute_segments(SPLICE_GRAPH_COMPLEX_DICT, TRANSCRIPTOME_COMPLEX_DICT))
        expected = SEGMENTS_COMPLEX
        self.assertEqual(actual, expected)



class TestComputeLinks(TestCase):
    """Tests for _compute_links"""

    def test_empty(self):
        """exfi.io.splice_graph_dict_to_gfa1._compute_links: empty case"""
        actual = list(_compute_links(SPLICE_GRAPH_EMPTY_DICT))
        expected = LINKS_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """exfi.io.splice_graph_dict_to_gfa1._compute_links: simple case"""
        actual = list(_compute_links(SPLICE_GRAPH_SIMPLE_DICT))
        expected = LINKS_SIMPLE
        self.assertEqual(actual, expected)

    def test_coplex(self):
        """exfi.io.splice_graph_dict_to_gfa1._compute_links: complex case"""
        actual = list(_compute_links(SPLICE_GRAPH_COMPLEX_DICT))
        expected = LINKS_COMPLEX
        self.assertEqual(actual, expected)



class TestComputeContainments(TestCase):
    """Tests for _compute_containments"""

    def test_empty(self):
        """exfi.io.splice_graph_dict_to_gfa1._compute_containments: empty case"""
        actual = list(_compute_containments(SPLICE_GRAPH_EMPTY_DICT))
        expected = CONTAINMENTS_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """exfi.io.splice_graph_dict_to_gfa1._compute_containments: simple case"""
        actual = list(_compute_containments(SPLICE_GRAPH_SIMPLE_DICT))
        expected = CONTAINMENTS_SIMPLE
        self.assertEqual(actual, expected)

    def test_coplex(self):
        """exfi.io.splice_graph_dict_to_gfa1._compute_containments: complex case"""
        actual = list(_compute_containments(SPLICE_GRAPH_COMPLEX_DICT))
        expected = CONTAINMENTS_COMPLEX
        self.assertEqual(actual, expected)



class TestComputePaths(TestCase):
    """Tests for exfi.io.splice_graph_dict_to_gfa1._compute_paths"""

    def test_empty(self):
        """exfi.io.splice_graph_dict_to_gfa1._compute_paths: empty case"""
        actual = list(_compute_paths(SPLICE_GRAPH_EMPTY_DICT))
        expected = PATHS_EMPTY
        self.assertEqual(actual, expected)

    def test_simple(self):
        """exfi.io.splice_graph_dict_to_gfa1._compute_paths: simple case"""
        actual = list(_compute_paths(SPLICE_GRAPH_SIMPLE_DICT))
        expected = PATHS_SIMPLE
        self.assertEqual(actual, expected)

    def test_complex(self):
        """exfi.io.splice_graph_dict_to_gfa1._compute_paths: complex case"""
        actual = list(_compute_paths(SPLICE_GRAPH_COMPLEX_DICT))
        expected = PATHS_COMPLEX
        self.assertEqual(actual, expected)


class TestSpliceGraphToGFA1(TestCase):
    """Tests for exfi.io.splice_graph_dict_to_gfa1.splice_graph_dict_to_gfa1"""

    def test_empty(self):
        """exfi.io.splice_graph_dict_to_gfa1.splice_graph_dict_to_gfa1: empty case"""
        tmp_file = tempfile.mkstemp()[1]
        splice_graph_dict_to_gfa1(
            splice_graph_dict=SPLICE_GRAPH_EMPTY_DICT,
            transcriptome_dict=TRANSCRIPTOME_EMPTY_DICT,
            filename=tmp_file
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            GFA_EMPTY_FN
        ))
        os.remove(tmp_file)

    def test_simple(self):
        """exfi.io.splice_graph_dict_to_gfa1.splice_graph_dict_to_gfa1: simple case"""
        tmp_file = tempfile.mkstemp()[1]
        splice_graph_dict_to_gfa1(
            splice_graph_dict=SPLICE_GRAPH_SIMPLE_DICT,
            transcriptome_dict=TRANSCRIPTOME_SIMPLE_DICT,
            filename=tmp_file #tmp_file
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            GFA_SIMPLE_FN
        ))
        os.remove(tmp_file)

    def test_multiple(self):
        """exfi.io.splice_graph_dict_to_gfa1.splice_graph_dict_to_gfa1: complex case"""
        tmp_file = tempfile.mkstemp()[1]
        splice_graph_dict_to_gfa1(
            splice_graph_dict=SPLICE_GRAPH_COMPLEX_DICT,
            transcriptome_dict=TRANSCRIPTOME_COMPLEX_DICT,
            filename=tmp_file
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            GFA_COMPLEX_FN
        ))
        os.remove(tmp_file)

if __name__ == "__main__":
    main()
