#!/usr/bin/env python3

"""
Tests for exfi.io.read_gfa1
"""

from unittest import TestCase, main

from exfi.io.read_gfa1 import \
        _overlap_str_to_int, \
        _process_segments, \
        _process_links, \
        _process_containments, \
        _process_paths, \
        read_gfa1


from tests.test_data import \
    GFA_EMPTY_FN, GFA_SIMPLE_FN, GFA_COMPLEX_FN, \
    SEGMENTS_EMPTY, SEGMENTS_SIMPLE, SEGMENTS_COMPLEX, \
    SEGMENTS_EMPTY_DICT, SEGMENTS_SIMPLE_DICT, SEGMENTS_COMPLEX_DICT, \
    LINKS_EMPTY, LINKS_SIMPLE, LINKS_COMPLEX, \
    LINKS_EMPTY_DICT, LINKS_SIMPLE_DICT, LINKS_COMPLEX_DICT, \
    CONTAINMENTS_EMPTY, CONTAINMENTS_SIMPLE, CONTAINMENTS_COMPLEX, \
    CONTAINMENTS_EMPTY_DICT, CONTAINMENTS_SIMPLE_DICT, CONTAINMENTS_COMPLEX_DICT, \
    PATHS_EMPTY, PATHS_SIMPLE, PATHS_COMPLEX, \
    PATHS_EMPTY_DICT, PATHS_SIMPLE_DICT, PATHS_COMPLEX_DICT



def _split_lines(list_of_elements):
    """Process line"""
    for element in list_of_elements:
        yield element.strip().split("\t")



class TestOverlapStrToInt(TestCase):
    """Tests of exfi.io.read_gfa1._overlap_str_to_int"""

    def test_match(self):
        """exfi.io.read_gfa1._overlap_str_to_int: wrong case"""
        self.assertEqual(_overlap_str_to_int("13M"), 13)

    def test_gap(self):
        """exfi.io.read_gfa1._overlap_str_to_int: wrong case"""
        self.assertEqual(_overlap_str_to_int("12G"), -12)

    def test_missing_letter(self):
        """exfi.io.read_gfa1._overlap_str_to_int: messy case"""
        with self.assertRaises(ValueError):
            _overlap_str_to_int("12")



class TestProcessSegments(TestCase):
    """Tests of exfi.io.read_gfa1._process_segments"""

    def test_empty(self):
        """exfi.io.read_gfa1._process_segments: empty case"""
        self.assertEqual(
            _process_segments(_split_lines(SEGMENTS_EMPTY)),
            SEGMENTS_EMPTY_DICT
        )

    def test_simple(self):
        """exfi.io.read_gfa1._process_segments: simple case"""
        self.assertEqual(
            _process_segments(_split_lines(SEGMENTS_SIMPLE)),
            SEGMENTS_SIMPLE_DICT
        )
    def test_complex(self):
        """exfi.io.read_gfa1._process_segments: complex case"""
        self.assertEqual(
            _process_segments(_split_lines(SEGMENTS_COMPLEX)),
            SEGMENTS_COMPLEX_DICT
        )



class TestProcessLinks(TestCase):
    """Tests of exfi.io.read_gfa1._process_links"""

    def test_empty(self):
        """exfi.io.read_gfa1._process_links: empty case"""
        self.assertEqual(
            _process_links(_split_lines(LINKS_EMPTY)),
            LINKS_EMPTY_DICT
        )

    def test_simple(self):
        """exfi.io.read_gfa1._process_links: simple case"""
        self.assertEqual(
            _process_links(_split_lines(LINKS_SIMPLE)),
            LINKS_SIMPLE_DICT
        )
    def test_complex(self):
        """exfi.io.read_gfa1._process_links: complex case"""
        self.assertEqual(
            _process_links(_split_lines(LINKS_COMPLEX)),
            LINKS_COMPLEX_DICT
        )



class TestProcessContainments(TestCase):
    """Tests of exfi.io.read_gfa1._process_containments"""

    def test_empty(self):
        """exfi.io.read_gfa1._process_containments: empty case"""
        self.assertEqual(
            _process_containments(_split_lines(CONTAINMENTS_EMPTY)),
            CONTAINMENTS_EMPTY_DICT
        )

    def test_simple(self):
        """exfi.io.read_gfa1._process_containments: simple case"""
        self.assertEqual(
            _process_containments(_split_lines(CONTAINMENTS_SIMPLE)),
            CONTAINMENTS_SIMPLE_DICT
        )
    def test_complex(self):
        """exfi.io.read_gfa1._process_containments: complex case"""
        self.assertEqual(
            _process_containments(_split_lines(CONTAINMENTS_COMPLEX)),
            CONTAINMENTS_COMPLEX_DICT
        )



class TestProcessPaths(TestCase):
    """Tests of exfi.io.read_gfa1._process_paths"""

    def test_empty(self):
        """exfi.io.read_gfa1._process_paths: empty case"""
        self.assertEqual(
            _process_paths(_split_lines(PATHS_EMPTY)),
            PATHS_EMPTY_DICT
        )

    def test_simple(self):
        """exfi.io.read_gfa1._process_paths: simple case"""
        self.assertEqual(
            _process_paths(_split_lines(PATHS_SIMPLE)),
            PATHS_SIMPLE_DICT
        )
    def test_complex(self):
        """exfi.io.read_gfa1._process_paths: complex case"""
        self.assertEqual(
            _process_paths(_split_lines(PATHS_COMPLEX)),
            PATHS_COMPLEX_DICT
        )


class TestReadGFA1(TestCase):
    """Tests for exfi.io.read_gfa1.read_gfa1"""

    def test_empty(self):
        """exfi.io.read_gfa1.read_gfa1: empty case"""
        self.assertEqual(
            read_gfa1(GFA_EMPTY_FN),
            {
                "header": "VN:Z:1.0",
                "segments": SEGMENTS_EMPTY_DICT,
                "links": LINKS_EMPTY_DICT,
                "containments": CONTAINMENTS_EMPTY_DICT,
                "paths": PATHS_EMPTY_DICT
            }
        )

    def test_simple(self):
        """exfi.io.read_gfa1.read_gfa1: simple case"""
        self.assertEqual(
            read_gfa1(GFA_SIMPLE_FN),
            {
                "header": "VN:Z:1.0",
                "segments": SEGMENTS_SIMPLE_DICT,
                "links": LINKS_SIMPLE_DICT,
                "containments": CONTAINMENTS_SIMPLE_DICT,
                "paths": PATHS_SIMPLE_DICT
            }
        )

    def test_complex(self):
        """exfi.io.read_gfa1.read_gfa1: complex case"""
        self.assertEqual(
            read_gfa1(GFA_COMPLEX_FN),
            {
                "header": "VN:Z:1.0",
                "segments": SEGMENTS_COMPLEX_DICT,
                "links": LINKS_COMPLEX_DICT,
                "containments": CONTAINMENTS_COMPLEX_DICT,
                "paths": PATHS_COMPLEX_DICT
            }
        )



if __name__ == "__main__":
    main()
