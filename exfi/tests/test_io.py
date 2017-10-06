#!/usr/bin/env python3

import unittest
from exfi.build_splicegraph import \
    build_splicegraph

from exfi.io import \
    write_gfa1

import networkx as nx
import pandas as pd
import numpy as np
from Bio import SeqIO
import filecmp
import tempfile
import os

index_simple = SeqIO.index(
    filename="exfi/tests/files/build_splicegraph/single.fa",
    format="fasta"
)

index_different = SeqIO.index(
    filename="exfi/tests/files/build_splicegraph/different_transcripts.fa",
    format="fasta"
)

path_simple = {"ENSDART00000161035.1": ["EXON00000000001"]}

path_different = {
    "ENSDART00000161035.1":
        ["EXON00000000001", "EXON00000000002", "EXON00000000003"]
    ,
    "ENSDART00000165342.1":
        ["EXON00000000004", "EXON00000000005", "EXON00000000006",
        "EXON00000000007", "EXON00000000008", "EXON00000000009",
        "EXON00000000010", "EXON00000000011", "EXON00000000012",
        "EXON00000000013", "EXON00000000014", "EXON00000000015"]
}

transcriptome_simple = SeqIO.index(
    filename="exfi/tests/files/build_splicegraph/transcriptome_simple.fa",
    format="fasta"
)

transcriptome_different = SeqIO.index(
    filename="exfi/tests/files/build_splicegraph/transcriptome_different.fa",
    format="fasta"
)


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
            "exfi/tests/files/build_splicegraph/empty.gfa"
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
            "exfi/tests/files/build_splicegraph/single.gfa"
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
            "exfi/tests/files/build_splicegraph/different.gfa"
        ))
        os.remove(tmp_file)


class TestGFA1ToGappedTranscript(unittest.TestCase):

    def test_empty(self):
        pass

    def test_simple(self):
        pass

    def test_multiple(self):
        pass



class TestGFA1ToExons(unittest.TestCase):

    def test_empty(self):
        pass

    def test_simple(self):
        pass

    def test_multiple(self):
        pass
