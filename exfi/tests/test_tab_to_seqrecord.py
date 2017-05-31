#!/usr/bin/env python3


from unittest import TestCase
from exfi.tab_to_seqrecord import tab_to_seqrecord
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class TestTabToSeqRecord(TestCase):

    def test_empty(self):
        """Test an empty iterable"""
        actual = [x for x in tab_to_seqrecord([])]
        expected = []
        self.assertEqual(actual, expected)

    def test_one(self):
        """Test one element"""
        case = ["tr1:1-10\tACACAGTCAGCTAGTCATCGTAGT"]
        actual = [x for x in tab_to_seqrecord(case)][0]
        expected = SeqRecord(
            id="tr1",
            seq="ACACAGTCAGCTAGTCATCGTAGT",
            description="tr1:1-10"
        )
        self.assertEqual(actual.id, expected.id)
        self.assertEqual(actual.seq, expected.seq)
        self.assertEqual(actual.description, expected.description)

    def test_multiple(self):
        """Test multiple elements"""
        case = [
            "tr1:1-10\tACACAGTCAGCTAGTCATCGTAGT",
            "tr1:15-20\tACACAGTCAGCTAGTCGATCGATCCATCGTAGT",
            "tr2:18-25\tACACAGTCAGCTAGTCATCG"
        ]
        actual = [x for x in tab_to_seqrecord(case)]
        expected = [
            SeqRecord(
                id="tr1",
                description="tr1:1-10",
                seq=Seq("ACACAGTCAGCTAGTCATCGTAGT")
            ),
            SeqRecord(
                id="tr1",
                description="tr1:15-20",
                seq=Seq("ACACAGTCAGCTAGTCGATCGATCCATCGTAGT")
            ),
            SeqRecord(
                id="tr2",
                description="tr2:18-25",
                seq=Seq("ACACAGTCAGCTAGTCATCG")
            )
        ]
        self.assertEqual(
            [x.id for x in actual],
            [x.id for x in expected]
        )
        self.assertEqual(
            [x.description for x in actual],
            [x.description for x in expected]
        )
        self.assertEqual(
            [x.seq for x in actual],
            [x.seq for x in expected]
        )
