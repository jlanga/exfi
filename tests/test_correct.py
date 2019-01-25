#!/usr/bin/env python3

"""
test_correct_splice_graph.py: tests for the exfi.correct submodule
"""



from unittest import TestCase, main

from tempfile import \
    mkstemp

from os import remove
from os.path import dirname

import pandas as pd

from exfi.build_baited_bloom_filter import \
    build_baited_bloom_filter

from exfi.find_exons import \
    find_exons

from exfi.io.fasta_to_dict import \
    fasta_to_dict

from exfi.io.bed import \
    BED4_COLS, BED4_DTYPES, \
    bed3_to_bed4

from exfi.correct import \
    prepare_sealer, \
    run_sealer, \
    collect_sealer_results, \
    apply_correction_to_bed4, \
    correct_bed4

from tests.custom_assertions import \
    CustomAssertions


def _compose_args(bloom_fn: str, gfa_fn: str) -> dict:
    """Compose a dict of args with two variables"""
    return {
        "kmer": 30,
        "bloom": bloom_fn,
        "bloom_size": "500M",
        "levels": 1,
        "fasta": "tests/correct/transcript.fa",
        "max_fp_bases": 5,
        "max_overlap": 10,
        "gfa1": gfa_fn,
        "threads": 4,
        "max_gap_size": 10,
        "reads": ["tests/correct/reads.fa"],
    }


TEMP = mkstemp()
TEMPDIR = dirname(TEMP[1])
TEMP_BLOOM = TEMP[1] + ".bloom"
TEMP_GFA = TEMP[1] + ".gfa"
ARGS = _compose_args(TEMP_BLOOM, TEMP_GFA)
build_baited_bloom_filter(ARGS)
BED3 = find_exons(ARGS)
BED4 = bed3_to_bed4(BED3)
TRANSCRIPTOME_DICT = fasta_to_dict(ARGS["fasta"])



SEALED_EDGES = pd.DataFrame(
    data=[
        ["ENSDART00000149335.2:485-1715", "ENSDART00000149335.2:1717-2286"],
        ["ENSDART00000149335.2:1717-2286", "ENSDART00000149335.2:2288-3379"]
    ],
    columns=["u", "v"]
)

BED4_CORRECTED = pd.DataFrame(
    data=[
        ["ENSDART00000149335.2", 0, 486, "ENSDART00000149335.2:0-486"],
        ["ENSDART00000149335.2", 485, 3379, "ENSDART00000149335.2:485-3379"]
    ],
    columns=BED4_COLS
).astype(BED4_DTYPES)


def tearDownModule():
    """Remove temporary bloom and temporary GFA files"""
    # pylint: disable=invalid-name
    remove(TEMP[1])
    remove(TEMP_BLOOM)



class TestPrepareSealer(TestCase, CustomAssertions):
    """_prepare_sealer(splice_graph, args)
    (nx.DiGraph, dict_of_parameters) -> str
    """

    def test_file_creation(self):
        """exfi.correct._prepare_sealer: test creation"""
        sealer_input_fn = prepare_sealer(BED4, TRANSCRIPTOME_DICT, ARGS)
        actual = fasta_to_dict(sealer_input_fn)
        expected = fasta_to_dict("tests/correct/to_seal.fa")
        remove(sealer_input_fn)
        self.assertEqualDict(actual, expected)



class TestRunSealer(TestCase, CustomAssertions):
    """run_sealer(sealer_input_fn, args):
    (str, dict) -> str
    """

    def test_run(self):
        """exfi.correct.run_sealer: test if runs"""
        sealer_in_fn = prepare_sealer(
            bed4=BED4, transcriptome_dict=TRANSCRIPTOME_DICT, args=ARGS
        )
        sealer_out_fn = run_sealer(sealer_input_fn=sealer_in_fn, args=ARGS)
        actual = fasta_to_dict("tests/correct/sealed.fa")
        expected = fasta_to_dict(sealer_out_fn)
        remove(sealer_in_fn)
        remove(sealer_out_fn)
        self.assertEqualDict(actual, expected)



class TestCollectSealerResults(TestCase):
    """collect_sealer_results(handle):
    (str) -> dict
    """

    def test_collect_empty(self):
        """exfi.correct.collect_sealer_results: empty case"""
        empty_file = mkstemp()
        sealer_output_fn = run_sealer(sealer_input_fn=empty_file[1], args=ARGS)
        observed = collect_sealer_results(filename=sealer_output_fn)
        remove(empty_file[1])
        remove(sealer_output_fn)
        print("shape = {}\n".format(observed.shape))
        self.assertTrue(observed.shape == (0, 2))

    def test_collect_somedata(self):
        """exfi.correct.collect_sealer_results: some data"""
        observed = collect_sealer_results(
            filename="tests/correct/sealed.fa"
        )
        print("observed:\n", observed)
        print("expected:\n", SEALED_EDGES)
        self.assertTrue(observed.equals(SEALED_EDGES))



class TestApplySealerCorrection(TestCase):
    """apply_correction_to_bed4(bed4, sealed_edges) -> bed4_corrected"""

    def test_empty_sealed(self):
        """exfi.correct.apply_correction_to_bed4: no sealing"""
        no_sealing = pd.DataFrame(columns=["u", "v"])
        observed = apply_correction_to_bed4(BED4, no_sealing)
        self.assertTrue(BED4.equals(observed))

    def test_some_data(self):
        """exfi.correct.apply_correction_to_bed4: no sealing"""
        observed = apply_correction_to_bed4(BED4, SEALED_EDGES)
        print("BED4:\n", BED4)
        print("Observed:\n", observed)
        print("Expected:\n", BED4_CORRECTED)
        self.assertTrue(observed.equals(BED4_CORRECTED))


class TestCorrectBED4(TestCase):
    """correct_bed4(bed4, transcriptome_dict, args) -> bed4_corrected"""

    def test_simple(self):
        """exfi.correct.correct_bed4: some data"""
        observed = correct_bed4(
            bed4=BED4, transcriptome_dict=TRANSCRIPTOME_DICT, args=ARGS
        )
        expected = BED4_CORRECTED
        self.assertTrue(observed.equals(expected))




if __name__ == '__main__':
    main()
    # Remove BF
