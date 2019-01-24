#!/usr/bin/env python3

"""tests.test_io.test_bed4_to_gfa1.py: tests for exfi.io.bed4_to_gfa1.py"""

from unittest import TestCase, main

from tempfile import mkstemp
import os
import filecmp

from exfi.io.bed4_to_gfa1 import \
    compute_header, \
    compute_segments, \
    compute_links, \
    compute_containments, \
    compute_paths, \
    bed4_to_gfa1

from tests.io.bed import \
    BED4_EMPTY, BED4_SIMPLE, BED4_COMPLEX

from tests.io.transcriptome_dicts import \
    TRANSCRIPTOME_EMPTY_DICT, TRANSCRIPTOME_SIMPLE_DICT, \
    TRANSCRIPTOME_COMPLEX_DICT

from tests.io.gfa1 import \
    HEADER, \
    SEGMENTS_EMPTY, SEGMENTS_SIMPLE, SEGMENTS_COMPLEX, \
    LINKS_EMPTY, LINKS_SIMPLE, LINKS_COMPLEX, \
    CONTAINMENTS_EMPTY, CONTAINMENTS_SIMPLE, CONTAINMENTS_COMPLEX, \
    PATHS_EMPTY, PATHS_SIMPLE, PATHS_COMPLEX, \
    GFA1_EMPTY_FN, GFA1_SIMPLE_FN, GFA1_COMPLEX_FN, \
    GFA1_COMPLEX_SOFT_FN, GFA1_COMPLEX_HARD_FN



class TestComputeHeader(TestCase):
    """Tests for exfi.io.bed4_to_gfa1.compute_header"""

    def test_header(self):
        """exfi.io.bed4_to_gfa1.compute_header: single test"""
        observed = compute_header()
        self.assertTrue(observed.equals(HEADER))



class TestComputeSegments(TestCase):
    """Tests for exfi.io.bed4_to_gfa1.compute_segments"""

    def test_empty(self):
        """exfi.io.bed4_to_gfa1.compute_segments: empty case"""
        observed = compute_segments(BED4_EMPTY, TRANSCRIPTOME_EMPTY_DICT)
        self.assertTrue(observed.equals(SEGMENTS_EMPTY))

    def test_simple(self):
        """exfi.io.bed4_to_gfa1.compute_segments: simple case"""
        observed = compute_segments(BED4_SIMPLE, TRANSCRIPTOME_SIMPLE_DICT)
        self.assertTrue(observed.equals(SEGMENTS_SIMPLE))

    def test_complex(self):
        """exfi.io.bed4_to_gfa1.compute_segments: complex case"""
        observed = compute_segments(BED4_COMPLEX, TRANSCRIPTOME_COMPLEX_DICT)
        self.assertTrue(observed.equals(SEGMENTS_COMPLEX))


class TestComputeLinks(TestCase):
    """Tests for exfi.io.bed4_to_gfa1.compute_links"""

    def test_empty(self):
        """exfi.io.bed4_to_gfa1.compute_links: empty case"""
        observed = compute_links(BED4_EMPTY)
        self.assertTrue(observed.equals(LINKS_EMPTY))

    def test_simple(self):
        """exfi.io.bed4_to_gfa1.compute_links: simple case"""
        observed = compute_links(BED4_SIMPLE)
        self.assertTrue(observed.equals(LINKS_SIMPLE))

    def test_complex(self):
        """exfi.io.bed4_to_gfa1.compute_links: complex case"""
        observed = compute_links(BED4_COMPLEX)
        # print("Observed", observed, observed.dtypes, sep="\n")
        # print("Expected", LINKS_COMPLEX,  LINKS_COMPLEX.dtypes, sep="\n")
        self.assertTrue(observed.equals(LINKS_COMPLEX))



class TestComputeContainments(TestCase):
    """Tests for exfi.io.bed4_to_gfa1.compute_containments"""

    def test_empty(self):
        """exfi.io.bed4_to_gfa1.compute_containments: empty case"""
        observed = compute_containments(BED4_EMPTY)
        self.assertTrue(observed.equals(CONTAINMENTS_EMPTY))

    def test_simple(self):
        """exfi.io.bed4_to_gfa1.compute_containments: simple case"""
        observed = compute_containments(BED4_SIMPLE)
        # print("Observed", observed, observed.dtypes, sep="\n")
        # print("Expected", CONTAINMENTS_SIMPLE,  CONTAINMENTS_SIMPLE.dtypes, sep="\n")
        self.assertTrue(observed.equals(CONTAINMENTS_SIMPLE))

    def test_complex(self):
        """exfi.io.bed4_to_gfa1.compute_containments: complex case"""
        observed = compute_containments(BED4_COMPLEX)
        # print("Observed", observed, observed.dtypes, sep="\n")
        # print("Expected", CONTAINMENTS_COMPLEX,  CONTAINMENTS_COMPLEX.dtypes, sep="\n")
        self.assertTrue(observed.equals(CONTAINMENTS_COMPLEX))



class TestComputePaths(TestCase):
    """Tests for exfi.io.bed4_to_gfa1.compute_paths"""

    def test_empty(self):
        """exfi.io.bed4_to_gfa1.compute_paths: empty case"""
        observed = compute_paths(BED4_EMPTY)
        print("Observed", observed, observed.dtypes, sep="\n")
        print("Expected", PATHS_EMPTY, PATHS_EMPTY.dtypes, sep="\n")
        self.assertTrue(observed.equals(PATHS_EMPTY))

    def test_simple(self):
        """exfi.io.bed4_to_gfa1.compute_paths: simple case"""
        observed = compute_paths(BED4_SIMPLE)
        print("Observed", observed, observed.dtypes, sep="\n")
        print("Expected", PATHS_SIMPLE, PATHS_SIMPLE.dtypes, sep="\n")
        self.assertTrue(observed.equals(PATHS_SIMPLE))

    def test_complex(self):
        """exfi.io.bed4_to_gfa1.compute_paths: complex case"""
        observed = compute_paths(BED4_COMPLEX)
        print("Observed", observed, observed.dtypes, sep="\n")
        print("Expected", PATHS_COMPLEX, PATHS_COMPLEX.dtypes, sep="\n")
        self.assertTrue(observed.equals(PATHS_COMPLEX))



class TestBED4TOGFA1(TestCase):
    """Tests for exfi.io.bed4_to_gfa1.bed4_to_gfa1"""

    def test_empty(self):
        """exfi.io.bed4_to_gfa1.bed4_to_gfa1: empty case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        bed4_to_gfa1(
            gfa1_fn=tmp_file,
            bed4=BED4_EMPTY,
            transcriptome_dict=TRANSCRIPTOME_EMPTY_DICT
        )
        self.assertTrue(filecmp.cmp(tmp_file, GFA1_EMPTY_FN))
        os.remove(tmp_file)

    def test_simple(self):
        """exfi.io.bed4_to_gfa1.bed4_to_gfa1: simple case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        bed4_to_gfa1(
            gfa1_fn=tmp_file,
            bed4=BED4_SIMPLE,
            transcriptome_dict=TRANSCRIPTOME_SIMPLE_DICT
        )
        self.assertTrue(filecmp.cmp(tmp_file, GFA1_SIMPLE_FN))
        os.remove(tmp_file)

    def test_complex(self):
        """exfi.io.bed4_to_gfa1.bed4_to_gfa1: complex case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        bed4_to_gfa1(
            gfa1_fn=tmp_file,
            bed4=BED4_COMPLEX,
            transcriptome_dict=TRANSCRIPTOME_COMPLEX_DICT
        )
        self.assertTrue(filecmp.cmp(tmp_file, GFA1_COMPLEX_FN))
        os.remove(tmp_file)

    def test_complex_soft(self):
        """exfi.io.bed4_to_gfa1.bed4_to_gfa1: complex soft masked case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        bed4_to_gfa1(
            gfa1_fn=tmp_file,
            bed4=BED4_COMPLEX,
            transcriptome_dict=TRANSCRIPTOME_COMPLEX_DICT,
            masking='soft'
        )
        self.assertTrue(filecmp.cmp(tmp_file, GFA1_COMPLEX_SOFT_FN))
        os.remove(tmp_file)

    def test_complex_hard(self):
        """exfi.io.bed4_to_gfa1.bed4_to_gfa1: complex hard masked case"""
        tmp_file = mkstemp()[1]
        print(tmp_file)
        bed4_to_gfa1(
            gfa1_fn=tmp_file,
            bed4=BED4_COMPLEX,
            transcriptome_dict=TRANSCRIPTOME_COMPLEX_DICT,
            masking='hard'
        )
        self.assertTrue(filecmp.cmp(tmp_file, GFA1_COMPLEX_HARD_FN))
        os.remove(tmp_file)




if __name__ == '__main__':
    main()
