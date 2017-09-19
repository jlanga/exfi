#!/usr/bin/env python3

from unittest import TestCase

from exfi.build_baited_bloom_filter import build_baited_bloom_filter

class TestBuildBaitedBloomFilter(TestCase):

    # def test_build_empty_transcriptome(self):
    #     '''build_baited_bloom_filter.py: try to bait with empty transcriptome
    #
    #     Note: biobloommaker fails
    #     '''
    #     build_baited_bloom_filter(
    #         transcriptome="exfi/tests/files/empty.txt",
    #         kmer=30,
    #         bloom_size="100M",
    #         levels=1,
    #         output_bloom="test.bf",
    #         threads=1,
    #         reads="exfi/bin/tests/files/empty.txt"
    #     )

    def test_build_empty_library(self):
        '''build_baited_bloom_filter.py: build a BF without reads'''
        build_baited_bloom_filter(
            transcriptome="exfi/tests/files/build_baited_bloom_filter/small_transcriptome.fa",
            kmer=30,
            bloom_size="100M",
            levels=1,
            output_bloom="test.bf",
            threads=1,
            reads=["exfi/bin/tests/files/empty.txt"]
        )

    def test_build_one_library(self):
        '''build_baited_bloom_filter.py: build the BF with one library'''
        build_baited_bloom_filter(
            transcriptome="exfi/tests/files/build_baited_bloom_filter/small_transcriptome.fa",
            kmer=30,
            bloom_size="100M",
            levels=1,
            output_bloom="test.bf",
            threads=1,
            reads=["exfi/bin/tests/build_baited_bloom_filter/reads_1.fq"]
        )

    def test_build_two_libraries(self):
        '''build_baited_bloom_filter.py: build the BF with two libraries'''
        build_baited_bloom_filter(
            transcriptome="exfi/tests/files/build_baited_bloom_filter/small_transcriptome.fa",
            kmer=30,
            bloom_size="100M",
            levels=1,
            output_bloom="test.bf",
            threads=1,
            reads=[
                "exfi/bin/tests/build_baited_bloom_filter/reads_1.fq",
                "exfi/bin/tests/build_baited_bloom_filter/reads_2.fq"
            ]
        )
