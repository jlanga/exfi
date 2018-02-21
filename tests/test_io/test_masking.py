#!/usr/bin/env python3

"""
Tests for exfi.io.masking
"""

from unittest import TestCase, main

from exfi.io.masking import \
    _process_overlap_cigar, \
    _soft_mask_right, \
    _soft_mask_left, \
    _soft_mask, \
    _hard_mask_right, \
    _hard_mask_left, \
    _hard_mask, \
    _mask

from tests.data import \
    OVERLAPS_COMPLEX

# pylint: disable=no-name-in-module
from tests.test_io.test_gfa1_to_exons import \
    EXONS_COMPLEX_DICT, EXONS_COMPLEX_SOFT_DICT, EXONS_COMPLEX_HARD_DICT

from tests.custom_assertions import CustomAssertions


class TestProcessOverlapCigar(TestCase):
    """Tests for exfi.io.masking._process_overlap_cigar"""

    def test_empty(self):
        """exfi.io.masking._process_overlap_cigar: empty case"""
        with self.assertRaises(IndexError):
            _process_overlap_cigar("")

    def test_wrong(self):
        """exfi.io.masking._process_overlap_cigar: wrong case"""
        with self.assertRaises(IndexError):
            _process_overlap_cigar("")

    def test_correct(self):
        """exfi.io.masking._process_overlap_cigar: correct case"""
        self.assertEqual(
            _process_overlap_cigar("13M"),
            ["M", 13]
        )



class TestSoftMaskRight(TestCase):
    """Tests for exfi.io.masking._soft_mask_right"""

    def test_soft_mask_right(self):
        """exfi.io.masking._soft_mask_right: simple"""
        actual = _soft_mask_right("AAAAA", 3)
        expected = "AAaaa"
        self.assertEqual(actual, expected)



class TestSoftMaskLeft(TestCase):
    """Tests for exfi.io.masking._soft_mask_left"""

    def test_soft_mask_left(self):
        """exfi.io.masking._soft_mask_left: simple"""
        actual = _soft_mask_left("AAAAA", 3)
        expected = "aaaAA"
        self.assertEqual(actual, expected)



class TestSoftMask(TestCase, CustomAssertions):
    """Tests for exfi.io.masking._soft_mask"""

    def test_soft_mask(self):
        """exfi.io.masking._soft_mask: simple"""
        actual = _soft_mask(EXONS_COMPLEX_DICT, OVERLAPS_COMPLEX)
        expected = EXONS_COMPLEX_SOFT_DICT
        self.assertEqualDict(actual, expected)



class TestHardMaskRight(TestCase):
    """Tests for exfi.io.masking._hard_mask_right"""

    def test_hard_mask_right(self):
        """exfi.io.masking._hard_mask_right: simple"""
        actual = _hard_mask_right("AAAAA", 3)
        expected = "AANNN"
        self.assertEqual(actual, expected)



class TestHardMaskLeft(TestCase):
    """Tests for exfi.io.masking._hard_mask_left"""

    def test_hard_mask_left(self):
        """exfi.io.masking._hard_mask_left: simple"""
        actual = _hard_mask_left("AAAAA", 3)
        expected = "NNNAA"
        self.assertEqual(actual, expected)



class TestHardMask(TestCase, CustomAssertions):
    """Tests for exfi.io.masking._hard_mask"""

    def test_hard_mask(self):
        """exfi.io.masking._hard_mask: simple"""
        actual = _hard_mask(EXONS_COMPLEX_DICT, OVERLAPS_COMPLEX)
        expected = EXONS_COMPLEX_HARD_DICT
        self.assertEqualDict(actual, expected)



class TestMask(TestCase, CustomAssertions):
    """Tests for exfi.io.masking._mask"""

    def test_no_mask(self):
        """exfi.io.masking._mask: no masking"""
        actual = _mask(EXONS_COMPLEX_DICT, OVERLAPS_COMPLEX, "none")
        expected = EXONS_COMPLEX_DICT
        self.assertEqualDict(actual, expected)

    def test_soft_mask(self):
        """exfi.io.masking._mask: soft masking"""
        actual = _mask(EXONS_COMPLEX_DICT, OVERLAPS_COMPLEX, "soft")
        expected = EXONS_COMPLEX_SOFT_DICT
        self.assertEqualDict(actual, expected)

    def test_hard_mask(self):
        """exfi.io.masking._mask: hard masking"""
        actual = _mask(EXONS_COMPLEX_DICT, OVERLAPS_COMPLEX, "hard")
        expected = EXONS_COMPLEX_HARD_DICT
        self.assertEqualDict(actual, expected)



if __name__ == "__main__":
    main()
