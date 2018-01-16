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



if __name__ == "__main__":
    main()
