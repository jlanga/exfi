#!/usr/bin/env python3

"""
Tests for exfi.io.masking
"""

from unittest import TestCase, main

from exfi.io import index_fasta

from exfi.io.masking import \
    _process_overlap_cigar, \
    _soft_mask_right, \
    _soft_mask_left, \
    _soft_mask, \
    _hard_mask_right, \
    _hard_mask_left, \
    _hard_mask, \
    _mask

from tests.test_data import \
    OVERLAPS_COMPLEX, \
    EXONS_COMPLEX_DICT, EXONS_COMPLEX_SOFT_DICT, EXONS_COMPLEX_HARD_DICT



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
        self.assertEqual(
            _soft_mask_right("AAAAA", 3),
            "AAaaa"
        )

class TestSoftMaskLeft(TestCase):
    """Tests for exfi.io.masking._soft_mask_left"""
    def test_soft_mask_left(self):
        """exfi.io.masking._soft_mask_left: simple"""
        self.assertEqual(
            _soft_mask_left("AAAAA", 3),
            "aaaAA"
        )

class TestSoftMask(TestCase):
    """Tests for exfi.io.masking._soft_mask"""
    def test_soft_mask(self):
        """exfi.io.masking._soft_mask: simple"""
        print(EXONS_COMPLEX_DICT)
        self.assertEqual(
            _soft_mask(EXONS_COMPLEX_DICT, OVERLAPS_COMPLEX),
            EXONS_COMPLEX_SOFT_DICT
        )


class TestHardMaskRight(TestCase):
    """Tests for exfi.io.masking._hard_mask_right"""
    def test_hard_mask_right(self):
        """exfi.io.masking._hard_mask_right: simple"""
        self.assertEqual(
            _hard_mask_right("AAAAA", 3),
            "AANNN"
        )

class TestHardMaskLeft(TestCase):
    """Tests for exfi.io.masking._hard_mask_left"""
    def test_hard_mask_left(self):
        """exfi.io.masking._hard_mask_left: simple"""
        self.assertEqual(
            _hard_mask_left("AAAAA", 3),
            "NNNAA"
        )

class TestHardMask(TestCase):
    """Tests for exfi.io.masking._hard_mask"""
    def test_hard_mask(self):
        """exfi.io.masking._hard_mask: simple"""
        self.assertEqual(
            _hard_mask(EXONS_COMPLEX_DICT, OVERLAPS_COMPLEX),
            EXONS_COMPLEX_HARD_DICT
        )

class TestMask(TestCase):
    """Tests for exfi.io.masking._mask"""
    def test_no_mask(self):
        """exfi.io.masking._mask: no masking"""
        self.assertEqual(
            _mask(EXONS_COMPLEX_DICT, OVERLAPS_COMPLEX, False, False),
            EXONS_COMPLEX_DICT
        )

    def test_soft_mask(self):
        """exfi.io.masking._mask: soft masking"""
        self.assertEqual(
            _mask(EXONS_COMPLEX_DICT, OVERLAPS_COMPLEX, True, False),
            EXONS_COMPLEX_SOFT_DICT
        )

    def test_hard_mask(self):
        """exfi.io.masking._mask: hard masking"""
        self.assertEqual(
            _mask(EXONS_COMPLEX_DICT, OVERLAPS_COMPLEX, False, True),
            EXONS_COMPLEX_HARD_DICT
        )



if __name__ == "__main__":
    main()
