#!/usr/bin/env python3

from unittest import TestCase
from exfi.build_baited_bloom_filter import \
    build_baited_bloom_filter
import tempfile
import os

class TestBuildBaitedBloomFilter(TestCase):

    @classmethod
    def test_build_empty_transcriptome(self):
        '''build_baited_bloom_filter.py: try to bait with empty transcriptome

        Note: biobloommaker fails
        '''
        folder = tempfile.mkdtemp()
        build_baited_bloom_filter(
            transcriptome="exfi/tests/files/empty.txt",
            kmer=30,
            bloom_size="100M",
            levels=1,
            output_bloom=folder + "test.bf",
            threads=1,
            reads="exfi/tests/files/empty.txt"
        )
        os.rmdir(folder)

    @classmethod
    def test_build_empty_library(self):
        '''build_baited_bloom_filter.py: build a BF without reads'''
        folder = tempfile.mkdtemp()
        build_baited_bloom_filter(
            transcriptome="exfi/tests/files/build_baited_bloom_filter/small_transcriptome.fa",
            kmer=30,
            bloom_size="100M",
            levels=1,
            output_bloom=folder + "test.bf",
            threads=1,
            reads=["exfi/tests/files/empty.txt"]
        )
        os.rmdir(folder)

    @classmethod
    def test_build_one_library(self):
        '''build_baited_bloom_filter.py: build the BF with one library'''
        folder = tempfile.mkdtemp()
        build_baited_bloom_filter(
            transcriptome="exfi/tests/files/build_baited_bloom_filter/small_transcriptome.fa",
            kmer=30,
            bloom_size="100M",
            levels=1,
            output_bloom=folder + "test.bf",
            threads=1,
            reads=["exfi/tests/files/build_baited_bloom_filter/reads_1.fq"]
        )
        os.rmdir(folder)

    @classmethod
    def test_build_two_libraries(self):
        '''build_baited_bloom_filter.py: build the BF with two libraries'''
        folder = tempfile.mkdtemp()
        build_baited_bloom_filter(
            transcriptome="exfi/tests/files/build_baited_bloom_filter/small_transcriptome.fa",
            kmer=30,
            bloom_size="100M",
            levels=1,
            output_bloom=folder + "test.bf",
            threads=1,
            reads=[
                "exfi/tests/files/build_baited_bloom_filter/reads_1.fq",
                "exfi/tests/files/build_baited_bloom_filter/reads_2.fq"
            ]
        )
        os.rmdir(folder)
