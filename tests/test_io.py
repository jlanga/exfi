#!/usr/bin/env python3

import unittest
from exfi.build_splicegraph import \
    build_splicegraph

from exfi.io import \
    write_gfa1, \
    gfa1_to_exons, \
    gfa1_to_gapped_transcript

import filecmp
import tempfile
import os

from tests.test_data import *


class TestWriteGFA1(unittest.TestCase):

    def test_empty(self):
        """Write an empty GFA1 (just header)"""
        tmp_file = tempfile.mkstemp()[1]
        write_gfa1(
            splice_graph=build_splicegraph({}),
            exons={},
            filename=tmp_file
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            empty_gfa
        ))
        os.remove(tmp_file)

    def test_simple(self):
        """Write a single exon GFA"""
        tmp_file = tempfile.mkstemp()[1]
        write_gfa1(
            splice_graph=build_splicegraph(index_simple),
            exons=index_simple,
            filename=tmp_file #tmp_file
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            single_gfa
        ))
        os.remove(tmp_file)

    def test_multiple(self):
        """Write a more complex GFA"""
        tmp_file = tempfile.mkstemp()[1]
        write_gfa1(
            splice_graph=build_splicegraph(index_different),
            exons=index_different,
            filename=tmp_file
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            different_gfa
        ))
        os.remove(tmp_file)



class TestGFA1ToExons(unittest.TestCase):

    def test_empty(self):
        """Convert an empty GFA1 to an empty exon FASTA"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=empty_gfa,
            fasta_out_fn=tmp_file,
            soft_mask_overlaps=False
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            empty_exons
        ))
        os.remove(tmp_file)


    def test_simple(self):
        """Convert a simple GFA1 to a single exon FASTA"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=single_gfa,
            fasta_out_fn=tmp_file,
            soft_mask_overlaps=False
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            single_exons
        ))
        os.remove(tmp_file)

    def test_multiple(self):
        """Convert a more complex GFA1 to multiple exon FASTA"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=different_gfa,
            fasta_out_fn=tmp_file,
            soft_mask_overlaps=False
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            different_exons
        ))
        os.remove(tmp_file)

    def test_multiple_soft(self):
        """Convert a more complex GFA1 to multiple soft masked exon FASTA"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=different_gfa,
            fasta_out_fn=tmp_file,
            soft_mask_overlaps=True
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            different_exons_soft
        ))
        os.remove(tmp_file)

    def test_multiple_hard(self):
        """Convert a more complex GFA1 to multiple hard masked exon FASTA"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=different_gfa,
            fasta_out_fn=tmp_file,
            hard_mask_overlaps=True
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            different_exons_hard
        ))
        os.remove(tmp_file)



class TestGFA1ToGappedTranscript(unittest.TestCase):

    def test_empty(self):
        """Convert an empty GFA1 to an empty gapped transcriptFASTA"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(
            gfa_in=empty_gfa,
            fasta_out=tmp_file
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            empty_gapped
        ))
        os.remove(tmp_file)

    def test_simple(self):
        """Convert a simple GFA1 to an one seq FASTA"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(
            gfa_in=single_gfa,
            fasta_out=tmp_file
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            single_gapped
        ))
        os.remove(tmp_file)

    def test_multiple(self):
        """Convert an more complex GFA1 to a mutli seq FASTA"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(
            gfa_in=different_gfa,
            fasta_out=tmp_file
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            different_gapped
        ))
        os.remove(tmp_file)

    def test_multiple_soft(self):
        """Convert an more complex GFA1 to a mutli seq FASTA with overlaps soft masked"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(
            gfa_in=different_gfa,
            fasta_out=tmp_file,
            soft_mask_overlaps=True
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            different_gapped_soft
        ))
        os.remove(tmp_file)


    def test_multiple_hard(self):
        """Convert an more complex GFA1 to a mutli seq FASTA with overlaps hard masked"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(
            gfa_in=different_gfa,
            fasta_out=tmp_file,
            number_of_ns=100,
            hard_mask_overlaps=True
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            different_gapped_hard
        ))
        os.remove(tmp_file)
