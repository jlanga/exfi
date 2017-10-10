#!/usr/bin/env python3

from unittest import TestCase
from exfi.build_baited_bloom_filter import \
    build_baited_bloom_filter
import tempfile
import shutil

def _simple_build_baited(transcriptome, reads, tmp_dir, tmp_bf):
    build_baited_bloom_filter(
        transcriptome=transcriptome,
        kmer=30,
        bloom_size="100M",
        levels=1,
        output_bloom=tmp_bf,
        threads=1,
        reads=reads
    )
    shutil.rmtree(tmp_dir)


class TestBuildBaitedBloomFilter(TestCase):

    @classmethod
    def test_build_empty_transcriptome(self):
        '''build_baited_bloom_filter.py: try to bait with empty transcriptome

        Note: biobloommaker fails
        '''
        transcriptome = "tests/files/empty.txt"
        reads = "tests/files/empty.txt"
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/empty_transcriptome.bf"
        _simple_build_baited(transcriptome, reads, tmp_dir, tmp_bf)

    @classmethod
    def test_build_empty_library(self):
        '''build_baited_bloom_filter.py: build a BF without reads'''
        transcriptome = "tests/files/build_baited_bloom_filter/small_transcriptome.fa"
        reads = ["tests/files/empty.txt"]
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/empty_library.bf"
        _simple_build_baited(transcriptome, reads, tmp_dir, tmp_bf)

    @classmethod
    def test_build_one_library(self):
        '''build_baited_bloom_filter.py: build the BF with one library'''
        transcriptome = "tests/files/build_baited_bloom_filter/small_transcriptome.fa"
        reads = ["tests/files/build_baited_bloom_filter/reads_1.fq"]
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/one_library.bf"
        _simple_build_baited(transcriptome, reads, tmp_dir, tmp_bf)

    @classmethod
    def test_build_two_libraries(self):
        '''build_baited_bloom_filter.py: build the BF with two libraries'''
        transcriptome = "tests/files/build_baited_bloom_filter/small_transcriptome.fa"
        reads = [
            "tests/files/build_baited_bloom_filter/reads_1.fq",
            "tests/files/build_baited_bloom_filter/reads_2.fq"
        ]
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/two_libraries.bf"
        _simple_build_baited(transcriptome, reads, tmp_dir, tmp_bf)
