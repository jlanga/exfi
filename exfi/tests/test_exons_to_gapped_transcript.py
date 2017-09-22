#!/usr/bin/env python3

from unittest import TestCase
from exfi.exons_to_gapped_transcript import \
    build_transcript_to_exon_dict, \
    exon_dict_to_gapped_transcript
from Bio import SeqIO
from exfi.tests.auxiliary_functions import \
    CustomAssertions, \
    _fasta_to_list

def fasta_to_transcript_to_exon_dict(fasta_fn):
    """Pipe to read fasta and build dict"""
    exons = _fasta_to_list(fasta_fn)
    return build_transcript_to_exon_dict(exons)

class TestBuildTranscriptToExonDict(TestCase):

    def test_empty_exons(self):
        actual = build_transcript_to_exon_dict([])
        expected = {}
        self.assertEqual(actual, expected)

    def test_single_exon(self):
        actual = fasta_to_transcript_to_exon_dict(
            "exfi/tests/files/exons_to_gapped_transcript/single.fa",
        )
        expected = {"ENSDART00000161035.1": ("EXON00000000001", )}  # TODO
        self.assertEqual(actual, expected)

    def test_ordered_exons(self):
        actual = fasta_to_transcript_to_exon_dict(
            "exfi/tests/files/exons_to_gapped_transcript/ordered.fa",
        )
        expected = {"ENSDART00000161035.1":
            ("EXON00000000001", "EXON00000000002", "EXON00000000003")
        }
        self.assertEqual(actual, expected)

    def test_disordered_exons(self):
        actual = fasta_to_transcript_to_exon_dict(
            "exfi/tests/files/exons_to_gapped_transcript/disordered.fa",
        )
        expected = {"ENSDART00000161035.1":
            ("EXON00000000001", "EXON00000000002", "EXON00000000003")
        }
        self.assertEqual(actual, expected)

    def test_different_transcripts(self):
        actual = fasta_to_transcript_to_exon_dict(
            "exfi/tests/files/exons_to_gapped_transcript/different_transcripts.fa",
        )
        expected = {
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
        self.assertEqual(actual, expected)


class TestExonsToGappedTranscript(TestCase, CustomAssertions):
    # exon_dict_to_gapped_transcript(transcript_to_exons, exome_fn, number_of_ns=100)
    def test_empty_both(self):
        transcript_to_exons = {}
        exome_fn = "exfi/tests/files/exons_to_gapped_transcript/empty_exome.fa"
        actual = list(exon_dict_to_gapped_transcript(
            transcript_to_exons,
            exome_fn
        ))
        expected = []
        self.assertEqual(actual, expected)

    def test_empty_exome(self):
        transcript_to_exons = {"ENSDART00000161035.1":
            ("EXON00000000001", "EXON00000000002", "EXON00000000003")
        }
        exome_fn = "exfi/tests/files/exons_to_gapped_transcript/empty_exome.fa"
        actual = list(exon_dict_to_gapped_transcript(
            transcript_to_exons,
            exome_fn
        ))
        expected = list(SeqIO.parse(
            handle="exfi/tests/files/exons_to_gapped_transcript/empty_exome_result.fa",
            format="fasta"
        ))
        self.assertEqualListOfSeqrecords(actual, expected)

    def test_empty_transcriptome(self):
        transcript_to_exons = {}
        exome_fn = "exfi/tests/files/exons_to_gapped_transcript/ordered.fa"
        actual = list(exon_dict_to_gapped_transcript(
            transcript_to_exons,
            exome_fn
        ))
        expected = []
        self.assertEqualListOfSeqrecords(actual, expected)

    def test_single_exon(self):
        transcript_to_exons = {"ENSDART00000161035.1" : ("EXON00000000001", )}
        exome_fn = "exfi/tests/files/exons_to_gapped_transcript/single.fa"
        actual = list(exon_dict_to_gapped_transcript(
            transcript_to_exons,
            exome_fn
        ))
        print(actual)
        expected = list(SeqIO.parse(
            handle="exfi/tests/files/exons_to_gapped_transcript/single_exon_result.fa",
            format="fasta"
        ))
        self.assertEqualListOfSeqrecords(actual, expected)

    def test_multiple(self):
        transcript_to_exons = {"ENSDART00000161035.1":
            ("EXON00000000001", "EXON00000000002", "EXON00000000003")
        }
        exome_fn = "exfi/tests/files/exons_to_gapped_transcript/ordered.fa"
        actual = list(exon_dict_to_gapped_transcript(
            transcript_to_exons,
            exome_fn
        ))
        expected = list(SeqIO.parse(
            handle="exfi/tests/files/exons_to_gapped_transcript/paste_ordered.fa",
            format="fasta"
        ))
        self.assertEqualListOfSeqrecords(actual, expected)
