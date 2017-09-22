#!/usr/bin/env python3

from unittest import TestCase
from exfi.exons_to_gapped_transcript import \
    build_transcript_to_exon_dict, \
    exon_dict_to_gapped_transcript
from Bio import SeqIO
from exfi.tests.auxiliary_functions import \
    CustomAssertions, \
    _fasta_to_list

def _fasta_to_transcript_to_exon_dict(fasta_fn):
    """Pipe to read fasta and build dict"""
    exons = _fasta_to_list(fasta_fn)
    return build_transcript_to_exon_dict(exons)

very_simple_dict = {"ENSDART00000161035.1": ("EXON00000000001", )}  # TODO
simple_dict = {
    "ENSDART00000161035.1": (
        "EXON00000000001", "EXON00000000002", "EXON00000000003"
    )
}
harder_dict ={
    "ENSDART00000161035.1": (
        "EXON00000000001", "EXON00000000002", "EXON00000000003"
    ),
    "ENSDART00000165342.1": (
        "EXON00000000004", "EXON00000000005", "EXON00000000006",
        "EXON00000000007", "EXON00000000008", "EXON00000000009",
        "EXON00000000010", "EXON00000000011", "EXON00000000012",
        "EXON00000000013", "EXON00000000014", "EXON00000000015"
    )
}


class TestBuildTranscriptToExonDict(TestCase):

    def test_empty_exons(self):
        """exons_to_gapped_transcript.py: build the transcript to exon dict
        given an empty input"""
        self.assertEqual(build_transcript_to_exon_dict([]), {})

    def test_single_exon(self):
        """exons_to_gapped_transcript.py: build the transcript to exon dict
        given a single exon"""
        self.assertEqual(
            _fasta_to_transcript_to_exon_dict(
                "exfi/tests/files/exons_to_gapped_transcript/single.fa",
            ),
            very_simple_dict  # TODO
        )

    def test_ordered_exons(self):
        """exons_to_gapped_transcript.py: build the transcript to exon dict
        given already ordered exons"""
        self.assertEqual(
            _fasta_to_transcript_to_exon_dict(
                "exfi/tests/files/exons_to_gapped_transcript/ordered.fa",
            ),
            simple_dict
        )

    def test_disordered_exons(self):
        """exons_to_gapped_transcript.py: build the transcript to exon dict
        given disordered exons"""
        self.assertEqual(
            _fasta_to_transcript_to_exon_dict(
                "exfi/tests/files/exons_to_gapped_transcript/disordered.fa",
            ),
            simple_dict
        )

    def test_different_transcripts(self):
        """exons_to_gapped_transcript.py: build the transcript to exon dict
        given exons from different transcripts"""
        self.assertEqual(
            _fasta_to_transcript_to_exon_dict(
                "exfi/tests/files/exons_to_gapped_transcript/different_transcripts.fa",
            ),
            harder_dict
        )


class TestExonsToGappedTranscript(TestCase, CustomAssertions):
    # exon_dict_to_gapped_transcript(transcript_to_exons, exome_fn, number_of_ns=100)
    def test_empty_both(self):
        exome_fn = "exfi/tests/files/exons_to_gapped_transcript/empty_exome.fa"
        self.assertEqual(
            list(exon_dict_to_gapped_transcript(
                {},
                exome_fn
            )),
            []
        )

    def test_empty_exome(self):
        exome_fn = "exfi/tests/files/exons_to_gapped_transcript/empty_exome.fa"
        actual = list(exon_dict_to_gapped_transcript(simple_dict, exome_fn))
        expected = _fasta_to_list(
            "exfi/tests/files/exons_to_gapped_transcript/empty_exome_result.fa"
        )
        self.assertEqualListOfSeqrecords(actual, expected)

    def test_empty_transcriptome(self):
        exome_fn = "exfi/tests/files/exons_to_gapped_transcript/ordered.fa"
        actual = list(exon_dict_to_gapped_transcript({}, exome_fn))
        self.assertEqualListOfSeqrecords(actual, [])

    def test_single_exon(self):
        exome_fn = "exfi/tests/files/exons_to_gapped_transcript/single.fa"
        actual = list(exon_dict_to_gapped_transcript(
            very_simple_dict, exome_fn
        ))
        expected = _fasta_to_list(
            "exfi/tests/files/exons_to_gapped_transcript/single_exon_result.fa"
        )
        self.assertEqualListOfSeqrecords(actual, expected)

    def test_multiple(self):
        exome_fn = "exfi/tests/files/exons_to_gapped_transcript/ordered.fa"
        actual = list(exon_dict_to_gapped_transcript(
            simple_dict, exome_fn
        ))
        expected = _fasta_to_list(
            "exfi/tests/files/exons_to_gapped_transcript/paste_ordered.fa",
        )
        self.assertEqualListOfSeqrecords(actual, expected)
