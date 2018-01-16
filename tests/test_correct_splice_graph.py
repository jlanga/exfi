#!/usr/bin/env python3

"""
test_correct_splice_graph.py: tests for the exfi.correct_splice_graph submodule
"""


from unittest import TestCase, main

from tempfile import \
    mkstemp

from os import remove

import networkx as nx

from Bio import \
    SeqIO

from exfi.build_baited_bloom_filter import \
    build_baited_bloom_filter

from exfi.find_exons import \
    _find_exons_pipeline

from exfi.build_splice_graph import \
    build_splice_graph

from exfi.correct_splice_graph import \
    _prepare_sealer, \
    _run_sealer, \
    _collect_sealer_results, \
    _sculpt_graph, \
    correct_splice_graph



from tests.auxiliary_functions import \
    CustomAssertions



TEMP_BLOOM = mkstemp()
TEMP_GFA = mkstemp()


ARGS = {
    "kmer": 30,
    "input_bloom": TEMP_BLOOM[1],
    "bloom_size": "500M",
    "levels": 1,
    "input_fasta": "tests/correct_splice_graph/transcript.fa",
    "max_fp_bases": 5,
    "max_overlap": 10,
    "output_gfa": TEMP_GFA[1],
    "threads": 4,
    "max_gap_size": 10,
    "reads": ["tests/correct_splice_graph/reads.fa"],
    "output_bloom" : TEMP_BLOOM[1],
    "bloom_filter": TEMP_BLOOM[1]
}

build_baited_bloom_filter(ARGS)

# Get predicted exons in bed format
POSITIVE_EXONS_BED = list(_find_exons_pipeline(ARGS))

# Build splice graph
SPLICE_GRAPH = build_splice_graph(POSITIVE_EXONS_BED)




class TestPrepareSealer(TestCase, CustomAssertions):
    """_prepare_sealer(splice_graph, args)
    (nx.DiGraph, dict_of_parameters) -> str
    """
    def test_file_creation(self):
        """exfi.correct_splice_graph._prepare_sealer: test creation"""
        sealer_input_fn = _prepare_sealer(SPLICE_GRAPH, ARGS)
        actual = list(SeqIO.parse(sealer_input_fn, format="fasta"))
        expected = list(SeqIO.parse("tests/correct_splice_graph/to_seal.fa", format="fasta"))
        self.assertEqualListOfSeqrecords(actual, expected)
        remove(sealer_input_fn)



class TestRunSealer(TestCase, CustomAssertions):
    """_run_sealer(sealer_input_fn, args):
    (str, dict) -> str
    """
    def test_run(self):
        """exfi.correct_splice_graph._run_sealer: test if runs"""
        sealer_in_fn = _prepare_sealer(SPLICE_GRAPH, ARGS)
        sealer_out_fn = _run_sealer(sealer_input_fn=sealer_in_fn, args=ARGS)
        self.assertEqualListOfSeqrecords(
            list(SeqIO.parse("tests/correct_splice_graph/sealed.fa", format="fasta")),
            list(SeqIO.parse(sealer_out_fn, "fasta"))
        )
        remove(sealer_in_fn)
        remove(sealer_out_fn)



class TestCollectSealerResults(TestCase):
    """_collect_sealer_results(handle):
    (str) -> dict
    """
    def test_collect_empty(self):
        """exfi.correct_splice_graph._collect_sealer_results: empty case"""
        empty_file = mkstemp()
        sealer_output_fn = _run_sealer(sealer_input_fn=empty_file[1], args=ARGS)
        # Collect sealer results
        edge2fill = _collect_sealer_results(handle=sealer_output_fn)
        self.assertEqual(edge2fill, {})

    def test_collect_somedata(self):
        """exfi.correct_splice_graph._collect_sealer_results: some data"""
        edge2fill = _collect_sealer_results(
            handle="tests/correct_splice_graph/sealed.fa"
        )
        self.assertEqual(
            edge2fill,
            {
                'ENSDART00000149335.2:1717-2286': 'ENSDART00000149335.2:2288-3379',
                'ENSDART00000149335.2:485-1715': 'ENSDART00000149335.2:1717-2286'
            }
        )


class TestSculptGraph(TestCase):
    """_sculpt_graph(splice_graph, edge2fill):
    (nx.DiGraph, dict) -> nx.DiGraph
    """

    def test_sculpt_empty_data(self):
        """exfi.correct_splice_graph._sculpt_graph: empty case"""
        sealed_graph = _sculpt_graph(SPLICE_GRAPH, {})
        self.assertTrue(nx.is_isomorphic(
            sealed_graph,
            SPLICE_GRAPH
        ))


    def test_scuplt_real_data(self):
        """exfi.correct_splice_graph._sculpt_graph: some data"""
        test_graph = nx.DiGraph()
        test_graph.add_edge(
            u="ENSDART00000149335.2:0-486",
            v="ENSDART00000149335.2:1717-2286"
        )
        edge2fill = _collect_sealer_results(
            handle="tests/correct_splice_graph/sealed.fa"
        )
        sealed_graph = _sculpt_graph(SPLICE_GRAPH, edge2fill)
        self.assertTrue(nx.is_isomorphic(
            sealed_graph,
            test_graph
        ))



class TestCorrectSpliceGraph(TestCase):
    """ correct_splice_graph(splice_graph, args):
    (nx.DiGraph, int) -> nx.DiGraph
    """
    def test_correct_splice_graph(self):
        """exfi.correct_splice_graph.correct_splice_graph: some data"""
        test_graph = nx.DiGraph()
        test_graph.add_edge(
            u="ENSDART00000149335.2:0-486",
            v="ENSDART00000149335.2:485-3379"
        )
        splice_graph = build_splice_graph(POSITIVE_EXONS_BED)
        sealed_graph = correct_splice_graph(splice_graph, ARGS)
        print(test_graph.nodes())
        print(test_graph.edges())
        print(sealed_graph.nodes())
        print(sealed_graph.edges())
        self.assertTrue(nx.is_isomorphic(
            sealed_graph,
            test_graph
        ))


if __name__ == '__main__':
    main()
    # Remove BF
    remove(TEMP_BLOOM, TEMP_GFA)
