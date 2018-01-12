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
        '''_simple_build_baited: try to bait with empty transcriptome

        Note: biobloommaker fails
        '''
        transcriptome = "tests/files/empty.txt"
        reads = "tests/files/empty.txt"
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/empty_transcriptome.bf"
        _simple_build_baited(transcriptome, reads, tmp_bf)
        self.assertTrue(isfile(tmp_bf))
        shutil.rmtree(tmp_dir)

    def test_build_empty_library(self):
        '''_simple_build_baited: build a BF without reads'''
        transcriptome = "tests/files/build_baited_bloom_filter/small_transcriptome.fa"
        reads = ["tests/files/empty.txt"]
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/empty_library.bf"
        _simple_build_baited(transcriptome, reads, tmp_bf)
        self.assertTrue(isfile(tmp_bf))
        shutil.rmtree(tmp_dir)

    def test_build_one_library(self):
        '''_simple_build_baited: build the BF with one library'''
        transcriptome = "tests/files/build_baited_bloom_filter/small_transcriptome.fa"
        reads = ["tests/files/build_baited_bloom_filter/reads_1.fq"]
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/one_library.bf"
        _simple_build_baited(transcriptome, reads, tmp_bf)
        self.assertTrue(isfile(tmp_bf))
        shutil.rmtree(tmp_dir)

    def test_build_two_libraries(self):
        '''_simple_build_baited: build the BF with two libraries'''
        transcriptome = "tests/files/build_baited_bloom_filter/small_transcriptome.fa"
        reads = [
            "tests/files/build_baited_bloom_filter/reads_1.fq",
            "tests/files/build_baited_bloom_filter/reads_2.fq"
        ]
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/two_libraries.bf"
        _simple_build_baited(transcriptome, reads, tmp_bf)
        self.assertTrue(isfile(tmp_bf))
        shutil.rmtree(tmp_dir)
        # shutil.rmtree(tmp_dir + "/categories_multiMatch.fa")
        # shutil.rmtree(tmp_dir + "/categories_noMatch.fa")
        # shutil.rmtree(tmp_dir + "/categories_summary.tsv")
        # shutil.rmtree(tmp_dir + "/categories_transcriptome.fa")
