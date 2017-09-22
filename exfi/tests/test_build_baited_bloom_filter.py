#!/usr/bin/env python3

from unittest import TestCase
from exfi.build_baited_bloom_filter import \
    build_baited_bloom_filter
import tempfile
import shutil

class TestBuildBaitedBloomFilter(TestCase):

    @classmethod
    def test_build_empty_transcriptome(self):
        '''build_baited_bloom_filter.py: try to bait with empty transcriptome

        Note: biobloommaker fails
        '''
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/test.bf"
        build_baited_bloom_filter(
            transcriptome="exfi/tests/files/empty.txt",
            kmer=30,
            bloom_size="100M",
            levels=1,
            output_bloom=tmp_bf,
            threads=1,
            reads="exfi/tests/files/empty.txt"
        )
        shutil.rmtree(tmp_dir)

    @classmethod
    def test_build_empty_library(self):
        '''build_baited_bloom_filter.py: build a BF without reads'''
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/test.bf"
        build_baited_bloom_filter(
            transcriptome="exfi/tests/files/build_baited_bloom_filter/small_transcriptome.fa",
            kmer=30,
            bloom_size="100M",
            levels=1,
            output_bloom=tmp_bf,
            threads=1,
            reads=["exfi/tests/files/empty.txt"]
        )
        shutil.rmtree(tmp_dir)

    @classmethod
    def test_build_one_library(self):
        '''build_baited_bloom_filter.py: build the BF with one library'''
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/test.bf"
        build_baited_bloom_filter(
            transcriptome="exfi/tests/files/build_baited_bloom_filter/small_transcriptome.fa",
            kmer=30,
            bloom_size="100M",
            levels=1,
            output_bloom=tmp_bf,
            threads=1,
            reads=["exfi/tests/files/build_baited_bloom_filter/reads_1.fq"]
        )
        shutil.rmtree(tmp_dir)

    @classmethod
    def test_build_two_libraries(self):
        '''build_baited_bloom_filter.py: build the BF with two libraries'''
        tmp_dir = tempfile.mkdtemp()
        tmp_bf = tmp_dir + "/test.bf"
        build_baited_bloom_filter(
            transcriptome="exfi/tests/files/build_baited_bloom_filter/small_transcriptome.fa",
            kmer=30,
            bloom_size="100M",
            levels=1,
            output_bloom=tmp_bf,
            threads=1,
            reads=[
                "exfi/tests/files/build_baited_bloom_filter/reads_1.fq",
                "exfi/tests/files/build_baited_bloom_filter/reads_2.fq"
            ]
        )
        shutil.rmtree(tmp_dir)
