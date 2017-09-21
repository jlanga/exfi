#!/usr/bin/env python3

from unittest import TestCase
from exfi.build_baited_bloom_filter import \
    build_baited_bloom_filter
import tempfile
import os
import time

to_delete = ["/tmp/" + report for report in [
    "categories_multiMatch.fa",
    "categories_noMatch.fa",
    "categories_summary.fa"
]]

# print(to_delete)

class TestBuildBaitedBloomFilter(TestCase):

    @classmethod
    def test_build_empty_transcriptome(self):
        '''build_baited_bloom_filter.py: try to bait with empty transcriptome

        Note: biobloommaker fails
        '''
        filename = tempfile.mktemp(suffix=".bf")
        build_baited_bloom_filter(
            transcriptome="exfi/tests/files/empty.txt",
            kmer=30,
            bloom_size="100M",
            levels=1,
            output_bloom=filename,
            threads=1,
            reads="exfi/tests/files/empty.txt"
        )
        # os.unlink(filename)  # no file is generated

    @classmethod
    def test_build_empty_library(self):
        '''build_baited_bloom_filter.py: build a BF without reads'''
        filename = tempfile.mktemp(suffix=".bf")
        build_baited_bloom_filter(
            transcriptome="exfi/tests/files/build_baited_bloom_filter/small_transcriptome.fa",
            kmer=30,
            bloom_size="100M",
            levels=1,
            output_bloom=filename,
            threads=1,
            reads=["exfi/tests/files/empty.txt"]
        )
        os.unlink(filename)

    @classmethod
    def test_build_one_library(self):
        '''build_baited_bloom_filter.py: build the BF with one library'''
        filename = tempfile.mktemp(suffix=".bf")
        build_baited_bloom_filter(
            transcriptome="exfi/tests/files/build_baited_bloom_filter/small_transcriptome.fa",
            kmer=30,
            bloom_size="100M",
            levels=1,
            output_bloom=filename,
            threads=1,
            reads=["exfi/tests/files/build_baited_bloom_filter/reads_1.fq"]
        )
        os.unlink(filename)

    @classmethod
    def test_build_two_libraries(self):
        '''build_baited_bloom_filter.py: build the BF with two libraries'''
        filename = tempfile.mktemp(suffix=".bf")
        build_baited_bloom_filter(
            transcriptome="exfi/tests/files/build_baited_bloom_filter/small_transcriptome.fa",
            kmer=30,
            bloom_size="100M",
            levels=1,
            output_bloom=filename,
            threads=1,
            reads=[
                "exfi/tests/files/build_baited_bloom_filter/reads_1.fq",
                "exfi/tests/files/build_baited_bloom_filter/reads_2.fq"
            ]
        )
        os.unlink(filename)
