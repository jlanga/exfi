#!/usr/bin/env python3

import unittest
from exfi.build_splicegraph import \
    build_splicegraph

from exfi.io import \
    write_gfa1, \
    gfa1_to_exons, \
    gfa1_to_gapped_transcript, \
    _clean_index

from Bio import SeqIO
import filecmp
import tempfile
import os



index_simple = _clean_index(SeqIO.index(
    filename="exfi/tests/files/build_splicegraph/single.fa",
    format="fasta"
))

index_different = _clean_index(SeqIO.index(
    filename="exfi/tests/files/build_splicegraph/different_transcripts.fa",
    format="fasta"
))

transcriptome_simple = _clean_index(SeqIO.index(
    filename="exfi/tests/files/build_splicegraph/transcriptome_simple.fa",
    format="fasta"
))

transcriptome_different = _clean_index(SeqIO.index(
    filename="exfi/tests/files/build_splicegraph/transcriptome_different.fa",
    format="fasta"
))


empty_gfa = "exfi/tests/files/io/empty.gfa"
single_gfa = "exfi/tests/files/io/single.gfa"
different_gfa = "exfi/tests/files/io/different.gfa"

empty_exons = "exfi/tests/files/io/empty_exons.fa"
single_exons = "exfi/tests/files/io/single_exons.fa"
different_exons = "exfi/tests/files/io/different_exons.fa"
different_exons_masked = "exfi/tests/files/io/different_exons_masked.fa"

empty_gapped = "exfi/tests/files/io/empty_gapped.fa"
single_gapped = "exfi/tests/files/io/single_gapped.fa"
different_gapped = "exfi/tests/files/io/different_gapped.fa"
different_gapped_masked = "exfi/tests/files/io/different_gapped_masked.fa"


class TestWriteGFA1(unittest.TestCase):

    def test_empty(self):
        """Write an empty GFA1 (just header)"""
        tmp_file = tempfile.mkstemp()[1]
        write_gfa1(
            splice_graph=build_splicegraph({}),
            transcript_index={},
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
            transcript_index=transcriptome_simple,
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
            splice_graph= build_splicegraph(index_different),
            transcript_index=transcriptome_different,
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
        """Convert an empty GFA1 to an empty FASTA"""
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

    def test_multiple_masked(self):
        """Convert a more complex GFA1 to multiple masked exon FASTA"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_exons(
            gfa_in_fn=different_gfa,
            fasta_out_fn=tmp_file,
            soft_mask_overlaps=True
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            different_exons_masked
        ))
        os.remove(tmp_file)

class TestGFA1ToGappedTranscript(unittest.TestCase):

    def test_empty(self):
        """Convert an empty GFA1 to an empty FASTA"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(
            gfa_in=empty_gfa,
            fasta_out=tmp_file,
            number_of_ns=100,
            soft_mask_overlaps=False
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
            fasta_out=tmp_file,
            number_of_ns=100,
            soft_mask_overlaps=False
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
            fasta_out=tmp_file,
            number_of_ns=100,
            soft_mask_overlaps=False
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            different_gapped
        ))
        os.remove(tmp_file)

    def test_multiple_masked(self):
        """Convert an more complex GFA1 to a mutli seq FASTA with overlaps masked"""
        tmp_file = tempfile.mkstemp()[1]
        gfa1_to_gapped_transcript(
            gfa_in=different_gfa,
            fasta_out=tmp_file,
            number_of_ns=100,
            soft_mask_overlaps=True
        )
        self.assertTrue(filecmp.cmp(
            tmp_file,
            different_gapped_masked
        ))
        os.remove(tmp_file)
