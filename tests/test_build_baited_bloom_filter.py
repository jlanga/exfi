#!/usr/bin/env python3

"""
Battery of tests for exfi.build_baited_bloom_filter
"""

from unittest import TestCase

import tempfile
import shutil
from os.path import isfile

from exfi.build_baited_bloom_filter import \
    build_baited_bloom_filter

READS_EMPTY_FN = "tests/build_baited_bloom_filter/reads_empty.fq"
READS_1_FN = "tests/build_baited_bloom_filter/reads_1.fq"
READS_2_FN = "tests/build_baited_bloom_filter/reads_2.fq"

TRANSCRIPTOME_EMPTY_FN = "tests/build_baited_bloom_filter/transcriptome_empty.fa"
TRANSCRIPTOME_SMALL_FN = "tests/build_baited_bloom_filter/transcriptome_small.fa"

def _simple_build_baited(transcriptome, reads, tmp_bf):
    """(str, list, str, str) -> None

    Simple wrapper for testing"""
    args = {
        "input_fasta": transcriptome,
        "kmer": 30,
        "bloom_size": "100M",
        "levels": 1,
        "output_bloom": tmp_bf,
        "threads": 1,
        "reads": reads
    }
    build_baited_bloom_filter(args)


class TestBuildBaitedBloomFilter(TestCase):
    """Tests for build_baited_bloom_filter functions"""

    def test_build_empty_transcriptome(self):
        '''exfi.build_baited_bloom_filter._simple_build_baited: try to bait with empty transcriptome

        Note: biobloommaker fails
        '''
        transcriptome = TRANSCRIPTOME_EMPTY_FN
        reads = READS_EMPTY_FN
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/empty_transcriptome.bf"
        _simple_build_baited(transcriptome, reads, tmp_bf)
        self.assertTrue(isfile(tmp_bf))
        shutil.rmtree(tmp_dir)

    def test_build_empty_library(self):
        '''exfi.build_baited_bloom_filter._simple_build_baited: build a BF without reads'''
        transcriptome = TRANSCRIPTOME_SMALL_FN
        reads = [READS_EMPTY_FN]
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/empty_library.bf"
        _simple_build_baited(transcriptome, reads, tmp_bf)
        self.assertTrue(isfile(tmp_bf))
        shutil.rmtree(tmp_dir)

    def test_build_one_library(self):
        '''exfi.build_baited_bloom_filter._simple_build_baited: build the BF with one library'''
        transcriptome = TRANSCRIPTOME_SMALL_FN
        reads = [READS_1_FN]
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/one_library.bf"
        _simple_build_baited(transcriptome, reads, tmp_bf)
        self.assertTrue(isfile(tmp_bf))
        shutil.rmtree(tmp_dir)

    def test_build_two_libraries(self):
        '''exfi.build_baited_bloom_filter._simple_build_baited: build the BF with two libraries'''
        transcriptome = TRANSCRIPTOME_SMALL_FN
        reads = [READS_1_FN, READS_2_FN]
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/two_libraries.bf"
        _simple_build_baited(transcriptome, reads, tmp_bf)
        self.assertTrue(isfile(tmp_bf))
        shutil.rmtree(tmp_dir)
