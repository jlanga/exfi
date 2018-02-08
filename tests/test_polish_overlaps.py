#!/usr/bin/env python3

"""tests.test_polish_overlaps.py: Tests for the exfi.polish_overlaps submodule"""


from unittest import \
    TestCase, \
    main

from exfi.polish_overlaps import \
    coord_add_left, \
    coord_add_right, \
    polish_overlaps, \
    polish_overlaps_dict

from tests.custom_assertions import \
    CustomAssertions

from tests.test_data import \
    SPLICE_GRAPH_EMPTY, SPLICE_GRAPH_SIMPLE, SPLICE_GRAPH_COMPLEX, \
    SPLICE_GRAPH_EMPTY_DICT, SPLICE_GRAPH_SIMPLE_DICT, SPLICE_GRAPH_COMPLEX_DICT, \
    TRANSCRIPTOME_EMPTY, TRANSCRIPTOME_SIMPLE, TRANSCRIPTOME_COMPLEX, \
    POLISHED_EMPTY, POLISHED_SIMPLE, POLISHED_COMPLEX, \
    POLISHED_EMPTY_DICT, POLISHED_SIMPLE_DICT, POLISHED_COMPLEX_DICT


ARGS = {
    "threads": 1
}

class TestCoordAddLeft(TestCase):
    """Tests for exfi.polish_overlaps.coord_add_left"""

    def test_empty(self):
        """exfi.polish_overlaps.coord_add_left: empty case"""
        with self.assertRaises(IndexError):
            coord_add_left([], 1)

    def test_short(self):
        """exfi.polish_overlaps.coord_add_left: short case"""
        with self.assertRaises(IndexError):
            coord_add_left([1], 1)

    def test_correct(self):
        """exfi.polish_overlaps.coord_add_left: correct case"""
        actual = coord_add_left(("tr", 1, 2), 2)
        expected = ("tr", 3, 2)
        self.assertEqual(actual, expected)

    def test_too_long(self):
        """exfi.polish_overlaps.coord_add_left: too long case"""
        actual = coord_add_left(("tr", 1, 2, 5, 5), 2)
        expected = ("tr", 3, 2, 5, 5)
        self.assertEqual(actual, expected)



class TestCoordAddRight(TestCase):
    """Tests for exfi.polish_overlaps.coord_add_right"""

    def test_empty(self):
        """exfi.polish_overlaps.coord_add_right: empty case"""
        with self.assertRaises(IndexError):
            coord_add_right([], 1)

    def test_short(self):
        """exfi.polish_overlaps.coord_add_right: short case"""
        with self.assertRaises(IndexError):
            coord_add_right([1], 1)

    def test_correct(self):
        """exfi.polish_overlaps.coord_add_right: correct case"""
        actual = coord_add_right(("tr", 1, 2), 2)
        expected = ("tr", 1, 4)
        self.assertEqual(actual, expected)

    def test_too_long(self):
        """exfi.polish_overlaps.coord_add_right: too long case"""
        actual = coord_add_right(("tr", 1, 2, 5, 5), 2)
        expected = ("tr", 1, 4, 5, 5)
        self.assertEqual(actual, expected)



class TestPolishOverlaps(TestCase, CustomAssertions):
    """Tests for exfi.polish_overlaps.polish_overlaps"""

    def test_empty(self):
        """exfi.polish_overlaps.polish_overlaps: empty case"""
        actual = polish_overlaps(SPLICE_GRAPH_EMPTY, TRANSCRIPTOME_EMPTY)
        expected = POLISHED_EMPTY
        self.assertEqualSpliceGraphs(actual, expected)

    def test_simple(self):
        """exfi.polish_overlaps.polish_overlaps: simple case"""
        actual = polish_overlaps(SPLICE_GRAPH_SIMPLE, TRANSCRIPTOME_SIMPLE)
        expected = POLISHED_SIMPLE
        self.assertEqualSpliceGraphs(actual, expected)

    def test_complex(self):
        """exfi.polish_overlaps.polish_overlaps: complex case"""
        actual = polish_overlaps(SPLICE_GRAPH_COMPLEX, TRANSCRIPTOME_COMPLEX)
        expected = POLISHED_COMPLEX
        self.assertEqualSpliceGraphs(actual, expected)



class TestPolishOverlapsDict(TestCase, CustomAssertions):
    """Tests for exfi.polish_overlaps.polish_overlaps_dict"""

    def test_empty(self):
        """exfi.polish_overlaps.polish_overlaps_dict: empty case"""
        actual = polish_overlaps_dict(SPLICE_GRAPH_EMPTY_DICT, TRANSCRIPTOME_EMPTY, ARGS)
        expected = POLISHED_EMPTY_DICT
        self.assertEqualDictOfSpliceGraphs(actual, expected)

    def test_simple(self):
        """exfi.polish_overlaps.polish_overlaps_dict: simple case"""
        actual = polish_overlaps_dict(SPLICE_GRAPH_SIMPLE_DICT, TRANSCRIPTOME_SIMPLE, ARGS)
        expected = POLISHED_SIMPLE_DICT
        self.assertEqualDictOfSpliceGraphs(actual, expected)

    def test_complex(self):
        """exfi.polish_overlaps.polish_overlaps_dict: complex case"""
        actual = polish_overlaps_dict(SPLICE_GRAPH_COMPLEX_DICT, TRANSCRIPTOME_COMPLEX, ARGS)
        expected = POLISHED_COMPLEX_DICT
        print(actual)
        self.assertEqualDictOfSpliceGraphs(actual, expected)



if __name__ == '__main__':
    main()
