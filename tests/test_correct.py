#!/usr/bin/env python3

"""
test_correct_splice_graph.py: tests for the exfi.correct submodule
"""


from unittest import TestCase, main

from tempfile import \
    mkstemp

from os import remove
from os.path import dirname

import networkx as nx

from exfi.build_baited_bloom_filter import \
    build_baited_bloom_filter

from exfi.find_exons import \
    _find_exons_pipeline

from exfi.build_splice_graph_dict import \
    build_splice_graph_dict

from exfi.io.fasta_to_dict import \
    fasta_to_dict

from exfi.correct import \
    _prepare_sealer, \
    _run_sealer, \
    _collect_sealer_results, \
    _filled_edges_by_transcript, \
    _rename_nodes_from_collapse, \
    _recompute_node2coord, \
    _recompute_edge2overlap, \
    _compute_new_node_ids, \
    _sculpt_graph, \
    correct_splice_graph_dict

from tests.custom_assertions import \
    CustomAssertions


def _compose_args(bloom_fn, gfa_fn):
    """Compose a dict of args with two variables"""
    return {
        "kmer": 30,
        "input_bloom": bloom_fn,
        "bloom_size": "500M",
        "levels": 1,
        "input_fasta": "tests/correct/transcript.fa",
        "max_fp_bases": 5,
        "max_overlap": 10,
        "output_gfa": gfa_fn,
        "threads": 4,
        "max_gap_size": 10,
        "reads": ["tests/correct/reads.fa"],
        "output_bloom" : bloom_fn,
        "bloom_filter": bloom_fn,
    }


TEMP = mkstemp()
TEMPDIR = dirname(TEMP[1])
TEMP_BLOOM = TEMP[1] + ".bloom"
TEMP_GFA = TEMP[1] + ".gfa"
ARGS = _compose_args(TEMP_BLOOM, TEMP_GFA)
build_baited_bloom_filter(ARGS)
POSITIVE_EXONS_BED = list(_find_exons_pipeline(ARGS))
SPLICE_GRAPH_DICT = build_splice_graph_dict(POSITIVE_EXONS_BED, ARGS)
SPLICE_GRAPH = SPLICE_GRAPH_DICT["ENSDART00000149335.2"]
EDGE2FILL = {
    ('ENSDART00000149335.2:485-1715', 'ENSDART00000149335.2:1717-2286'),
    ('ENSDART00000149335.2:1717-2286', 'ENSDART00000149335.2:2288-3379')
}
FILLED_EDGE_BY_TRANSCRIPT = {
    'ENSDART00000149335.2': EDGE2FILL
}


def partition(node_u, node_v, edge2fill):
    """Define partitions as how the graph should be filled"""
    graph = nx.DiGraph()
    graph.add_edges_from(edge2fill)
    if node_u in graph.nodes() and \
        node_v in graph.nodes() and \
        nx.has_path(G=graph, source=node_u, target=node_v):
        return True
    return False

FULL_PARTITION = lambda u, v: partition(u, v, EDGE2FILL)
COLLAPSED_GRAPH = nx.quotient_graph(SPLICE_GRAPH, partition=FULL_PARTITION)

QUOTIENT_RELABELING = {
    frozenset({
        'ENSDART00000149335.2:0-486'
    }):
        'ENSDART00000149335.2:0-486',
    frozenset({
        'ENSDART00000149335.2:485-1715',
        'ENSDART00000149335.2:1717-2286',
        'ENSDART00000149335.2:2288-3379'
    }): (
        'ENSDART00000149335.2:485-1715',
        'ENSDART00000149335.2:1717-2286',
        'ENSDART00000149335.2:2288-3379'
    )
}

QUOTIENT_RELABELED = nx.relabel_nodes(
    copy=True,
    G=COLLAPSED_GRAPH,
    mapping=QUOTIENT_RELABELING
)

NEW_NODE2COORD = {
    'ENSDART00000149335.2:0-486':
        (('ENSDART00000149335.2', 0, 486),),
    ('ENSDART00000149335.2:485-1715', 'ENSDART00000149335.2:1717-2286',
     'ENSDART00000149335.2:2288-3379'):
        (('ENSDART00000149335.2', 485, 3379),)
}

NEW_EDGE2OVERLAP = {
    ('ENSDART00000149335.2:485-1715',
     'ENSDART00000149335.2:1717-2286',
     'ENSDART00000149335.2:2288-3379'):
        (('ENSDART00000149335.2', 485, 3379),),
    'ENSDART00000149335.2:0-486':
        (('ENSDART00000149335.2', 0, 486),)
}

NEW_NODE_IDS = {
    ('ENSDART00000149335.2:485-1715', 'ENSDART00000149335.2:1717-2286',
     'ENSDART00000149335.2:2288-3379'):
        'ENSDART00000149335.2:485-3379',
    'ENSDART00000149335.2:0-486':
        'ENSDART00000149335.2:0-486'
}

SEALED_GRAPH = nx.DiGraph()
SEALED_GRAPH.add_nodes_from(["ENSDART00000149335.2:0-486", "ENSDART00000149335.2:485-3379"])
nx.set_node_attributes(SEALED_GRAPH, name="coordinates", values={
    "ENSDART00000149335.2:0-486": (("ENSDART00000149335.2", 0, 486),),
    "ENSDART00000149335.2:485-3379": (("ENSDART00000149335.2", 485, 3379),)
})
SEALED_GRAPH.add_edge(
    u="ENSDART00000149335.2:0-486",
    v="ENSDART00000149335.2:485-3379"
)
nx.set_edge_attributes(SEALED_GRAPH, name="overlaps", values={
    ("ENSDART00000149335.2:0-486", "ENSDART00000149335.2:485-3379"): 1
})

SEALED_GRAPH_DICT = {"ENSDART00000149335.2": SEALED_GRAPH}

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
        sealer_input_fn = _prepare_sealer(SPLICE_GRAPH_DICT, ARGS)
        actual = fasta_to_dict(sealer_input_fn)
        expected = fasta_to_dict("tests/correct/to_seal.fa")
        remove(sealer_input_fn)
        self.assertEqual(actual, expected)



class TestRunSealer(TestCase, CustomAssertions):
    """_run_sealer(sealer_input_fn, args):
    (str, dict) -> str
    """

    def test_run(self):
        """exfi.correct._run_sealer: test if runs"""
        sealer_in_fn = _prepare_sealer(SPLICE_GRAPH_DICT, ARGS)
        sealer_out_fn = _run_sealer(sealer_input_fn=sealer_in_fn, args=ARGS)
        actual = fasta_to_dict("tests/correct/sealed.fa")
        expected = fasta_to_dict(sealer_out_fn)
        remove(sealer_in_fn)
        remove(sealer_out_fn)
        self.assertEqual(actual, expected)



class TestCollectSealerResults(TestCase):
    """_collect_sealer_results(handle):
    (str) -> dict
    """

    def test_collect_empty(self):
        """exfi.correct._collect_sealer_results: empty case"""
        empty_file = mkstemp()
        sealer_output_fn = _run_sealer(sealer_input_fn=empty_file[1], args=ARGS)
        edge2fill = _collect_sealer_results(handle=sealer_output_fn)
        remove(empty_file[1])
        remove(sealer_output_fn)
        self.assertEqual(edge2fill, set())

    def test_collect_somedata(self):
        """exfi.correct._collect_sealer_results: some data"""
        edge2fill = _collect_sealer_results(
            handle="tests/correct/sealed.fa"
        )
        self.assertEqual(edge2fill, EDGE2FILL)



class TestFilledEdgeByTranscript(TestCase):
    """Tests for _filled_edges_by_transcript(splice_graph: nx.DiGraph, filled_edges: str) -> dict"""

    def test_empty(self):
        """exfi.correct._filled_edge_by_transcript: empty case"""
        initial = {}
        actual = _filled_edges_by_transcript(filled_edges=initial)
        expected = {}
        self.assertEqual(actual, expected)

    def test_some_data(self):
        """exfi.correct._filled_edge_by_transcript: some data"""
        initial = _collect_sealer_results(
            handle="tests/correct/sealed.fa"
        )
        actual = _filled_edges_by_transcript(filled_edges=initial)
        expected = FILLED_EDGE_BY_TRANSCRIPT
        self.assertEqual(actual, expected)



class TestRenameNodesFromCollapse(TestCase):
    """Tests for _rename_nodes_from_collapse

    _rename_nodes_from_collapse(quotient_graph: nx.DiGraph) -> dict
    """

    def test_empty(self):
        """exfi.correct._rename_nodes_from_collapse: empty case"""
        initial = nx.DiGraph()
        actual = _rename_nodes_from_collapse(initial)
        expected = {}
        self.assertEqual(actual, expected)

    def test_some_data(self):
        """exfi.correct._rename_nodes_from_collapse: some data"""
        actual = _rename_nodes_from_collapse(COLLAPSED_GRAPH)
        expected = QUOTIENT_RELABELING
        self.assertEqual(actual, expected)




class TestRecomputeNode2Coord(TestCase):
    """Tests for exfi.correct._recompute_node2coord

    _recompute_node2coord(component: nx.DiGraph, quotient_relabeled: nx.DiGraph) -> dict
    """

    def test_empty(self):
        """exfi.correct._recompute_node2coord: empty case"""
        actual = _recompute_node2coord(nx.DiGraph(), nx.DiGraph())
        expected = {}
        self.assertEqual(actual, expected)

    def test_some_data(self):
        """exfi.correct._recompute_node2coord: some data"""
        actual = _recompute_node2coord(
            component=SPLICE_GRAPH,
            quotient_relabeled=QUOTIENT_RELABELED
        )
        expected = NEW_NODE2COORD
        self.assertEqual(actual, expected)



class TestRecomputeEdge2Overlap(TestCase):
    """Tests for _recompute_edge2overlap

    _recompute_edge2overlap(component: nx.DiGraph, quotient_relabeled: nx.DiGraph) -> dict
    """

    def test_empty(self):
        """exfi.correct._recompute_edge2overlap: empty case"""
        actual = _recompute_edge2overlap(nx.DiGraph(), nx.DiGraph())
        expected = {}
        self.assertEqual(actual, expected)

    def test_some_data(self):
        """exfi.correct._recompute_edge2overlap: some data"""
        actual = _recompute_node2coord(
            component=SPLICE_GRAPH,
            quotient_relabeled=QUOTIENT_RELABELED
        )
        expected = NEW_EDGE2OVERLAP
        self.assertEqual(actual, expected)



class TestComputeNewNodeIds(TestCase):
    """Test for exfi.correct._compute_new_node_ids

    _compute_new_node_ids(quotient_relabeled: nx.DiGraph, component: nx.DiGraph) -> dict
    """
    def test_empty(self):
        """exfi.correct._compute_new_node_ids: empty case"""
        actual = _compute_new_node_ids(nx.DiGraph(), nx.DiGraph())
        expected = {}
        self.assertEqual(actual, expected)

    def test_some_data(self):
        """exfi.correct._compute_new_node_ids: some data"""
        actual = _compute_new_node_ids(QUOTIENT_RELABELED, SPLICE_GRAPH)
        expected = NEW_NODE_IDS
        self.assertEqual(actual, expected)



class TestSculptGraph(TestCase):
    """Tests for exfi.correct._sculpt_graph

    _sculpt_graph(splice_graph: nx.DiGraph, filled_edges: set) -> nx.DiGraph
    """

    def test_sculpt_empty_data(self):
        """exfi.correct._sculpt_graph: empty case"""
        sealed_graph = _sculpt_graph(SPLICE_GRAPH, {})
        self.assertTrue(nx.is_isomorphic(
            sealed_graph,
            SPLICE_GRAPH
        ))


    def test_scuplt_real_data(self):
        """exfi.correct._sculpt_graph: some data"""
        test_graph = nx.DiGraph()
        test_graph.add_edge(
            u="ENSDART00000149335.2:0-486",
            v="ENSDART00000149335.2:485-3379"
        )
        edge2fill = _collect_sealer_results(
            handle="tests/correct/sealed.fa"
        )
        sealed_graph = _sculpt_graph(
            splice_graph=SPLICE_GRAPH, filled_edges=edge2fill
        )
        self.assertTrue(nx.is_isomorphic(
            sealed_graph,
            test_graph
        ))



class TestCorrectSpliceGraphDict(TestCase, CustomAssertions):
    """Tests for exfi.correct_dict.correct_splice_graph_dict"""

    def test_correct_splice_graph_dict(self):
        """exfi.correct.correct_splice_graph: some data"""

        splice_graph_dict = build_splice_graph_dict(POSITIVE_EXONS_BED, ARGS)
        sealed_graph_dict = correct_splice_graph_dict(splice_graph_dict, ARGS)
        self.assertEqualDictOfSpliceGraphs(sealed_graph_dict, SEALED_GRAPH_DICT)



if __name__ == '__main__':
    main()
    # Remove BF
