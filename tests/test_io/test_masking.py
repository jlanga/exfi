#!/usr/bin/env python3

"""
Tests for exfi.io.masking
"""

from unittest import TestCase, main

from exfi.io.masking import \
    mask

from tests.io.bed import \
    EDGE2OVERLAP_COMPLEX, NODE2SEQUENCE_COMPLEX, NODE2SEQUENCE_COMPLEX_HARD, \
    NODE2SEQUENCE_COMPLEX_SOFT



class TestMask(TestCase):
    """Tests for exfi.io.masking._mask"""

    def test_no_mask(self):
        """exfi.io.masking._mask: no masking"""
        observed = mask(NODE2SEQUENCE_COMPLEX, EDGE2OVERLAP_COMPLEX, "none")
        print(
            "Observed:", observed, "",
            "Expected:", NODE2SEQUENCE_COMPLEX,
            sep="\n"
        )
        self.assertTrue(observed.equals(NODE2SEQUENCE_COMPLEX))

    def test_soft_mask(self):
        """exfi.io.masking._mask: soft masking"""
        observed = mask(NODE2SEQUENCE_COMPLEX, EDGE2OVERLAP_COMPLEX, "soft")
        print(
            "Observed:", observed, "",
            "Expected:", NODE2SEQUENCE_COMPLEX_SOFT,
            sep="\n"
        )
        self.assertTrue(observed.equals(NODE2SEQUENCE_COMPLEX_SOFT))

    def test_hard_mask(self):
        """exfi.io.masking._mask: hard masking"""
        observed = mask(NODE2SEQUENCE_COMPLEX, EDGE2OVERLAP_COMPLEX, "hard")
        print(
            "Observed:", observed, "",
            "Expected:", NODE2SEQUENCE_COMPLEX_HARD,
            sep="\n"
        )
        self.assertTrue(observed.equals(NODE2SEQUENCE_COMPLEX_HARD))



if __name__ == "__main__":
    main()
