from unittest import TestCase

from exfi.correct_splice_graph import \
    _coordinates_to_variables, \
    _prepare_sealer, \
    _run_sealer, \
    _collect_sealer_results , \
    _sculpt_graph, \
    correct_splice_graph

from exfi.build_baited_bloom_filter import \
    build_baited_bloom_filter

from exfi.find_exons import \
    _find_exons_pipeline, \
    _get_fasta

from exfi.build_splice_graph import \
    build_splice_graph

from Bio import \
    SeqIO

import networkx as nx

from tempfile import \
    mkstemp

from tests.auxiliary_functions import \
    CustomAssertions

from os import remove

temp_bloom = mkstemp()
temp_gfa = mkstemp()


args = {
    "kmer": 30,
    "bloom_filter": temp_bloom[1],
    "bloom_size": "500M",
    "levels": 1,
    "input_fasta": "tests/files/correct_splice_graph/transcript.fa",
    "max_fp_bases": 5,
    "max_overlap": 10,
    "output_gfa": temp_gfa[1],
    "threads": 4,
    "max_gap_size": 10,
    "reads": ["tests/files/correct_splice_graph/reads.fa"]
}

build_baited_bloom_filter(
    transcriptome=args["input_fasta"],
    kmer=args["kmer"],
    bloom_size=args["bloom_size"],
    levels=args["levels"],
    output_bloom=args["bloom_filter"],
    threads=args["threads"],
    reads=args["reads"]
)

# Get predicted exons in bed format
positive_exons_bed = _find_exons_pipeline(
    kmer=args["kmer"],
    bloom_filter_fn=args["bloom_filter"],
    transcriptome_fn=args["input_fasta"],
    max_fp_bases=args["max_fp_bases"],
    max_overlap=args["max_overlap"]
)

# Bed -> fasta
transcriptome_index = SeqIO.index(
    filename=args["input_fasta"],
    format="fasta"
)

# Build splice graph
splice_graph = build_splice_graph(positive_exons_bed)





class TestCoordinatesToVariables(TestCase):
    """ _coordinates_to_variables(coordinates):
    (string) -> (string, int, int)
    """
    def test_empty_string(self):
        """_coordinates_to_variables: case empty"""
        with self.assertRaises(ValueError):
            _coordinates_to_variables("")

    def test_incorrect_string(self):
        """_coordinates_to_variables: case incorrect"""
        with self.assertRaises(IndexError):
            _coordinates_to_variables("ENSDART00000161035:1")

    def test_correct_string(self):
        """_coordinates_to_variables: case correct"""
        self.assertEqual(
            _coordinates_to_variables("ENSDART00000161035:1-15"),
            ("ENSDART00000161035", 1, 15)
        )

    def test_messy_string(self):
        """_coordinates_to_variables: messy case"""
        self.assertEqual(
            _coordinates_to_variables("TRINITY_g14_c15_i5:1:2-15-12:1-15"),
            ("TRINITY_g14_c15_i5:1:2-15-12", 1, 15)
        )



class TestPrepareSealer(TestCase, CustomAssertions):
    """_prepare_sealer(splice_graph, args)
    (nx.DiGraph, dict_of_parameters) -> str
    """
    def test_file_creation(self):
        """_prepare_sealer: test creation"""
        print(splice_graph.nodes())
        print(splice_graph.edges())
        sealer_input_fn = _prepare_sealer(splice_graph, args)
        print(sealer_input_fn)
        actual = list(SeqIO.parse(
            sealer_input_fn,
            format="fasta"
        ))
        expected = list(SeqIO.parse(
            "tests/files/correct_splice_graph/to_seal.fa",
            format="fasta"
        ))
        print(actual)
        print(expected)
        self.assertEqualListOfSeqrecords(
            actual,
            expected

        )
        remove(sealer_input_fn)



class TestRunSealer(TestCase, CustomAssertions):
    """_run_sealer(sealer_input_fn, args):
    (str, dict) -> str
    """
    def test_run(self):
        """_run_sealer: test if runs"""
        sealer_in_fn = _prepare_sealer(splice_graph, args)
        sealer_out_fn = _run_sealer(
            sealer_input_fn=sealer_in_fn,
            args=args
        )
        self.assertEqualListOfSeqrecords(
            list(SeqIO.parse(
                "tests/files/correct_splice_graph/sealed.fa",
                format="fasta"
            )),
            list(SeqIO.parse(
                sealer_out_fn, "fasta"
            ))
        )
        remove(sealer_in_fn)
        remove(sealer_out_fn)



class TestCollectSealerResults(TestCase):
    """_collect_sealer_results(handle):
    (str) -> dict
    """
    def test_collect_empty(self):
        """_collect_sealer_results: empty case"""
        empty_file = mkstemp()
        sealer_output_fn = _run_sealer(
            sealer_input_fn=empty_file[1],
            args=args
        )
        # Collect sealer results
        edge2fill = _collect_sealer_results(handle=sealer_output_fn)
        self.assertEqual(edge2fill, {})

    def test_collect_somedata(self):
        """_collect_sealer_results: some data"""
        edge2fill = _collect_sealer_results(
            handle="tests/files/correct_splice_graph/sealed.fa"
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
        """_sculpt_graph: empty case"""
        sealed_graph = _sculpt_graph(splice_graph, {}, transcriptome_index)
        self.assertTrue(nx.is_isomorphic(
            sealed_graph,
            splice_graph
        ))


    def test_scuplt_real_data(self):
        """_sculpt_graph: some data"""
        test_graph = nx.DiGraph()
        test_graph.add_edge(
            u="ENSDART00000149335.2:0-486",
            v="ENSDART00000149335.2:1717-2286"
        )
        edge2fill = _collect_sealer_results(
            handle="tests/files/correct_splice_graph/sealed.fa"
        )
        sealed_graph = _sculpt_graph(splice_graph, edge2fill, transcriptome_index)
        self.assertTrue(nx.is_isomorphic(
            sealed_graph,
            test_graph
        ))



class TestCorrectsplice_graph(TestCase):
    """ correct_splice_graph(splice_graph, args):
    (nx.DiGraph, int) -> nx.DiGraph
    """
    def test_correct_splice_graph(self):
        """correct_splice_graph: some data"""
        test_graph = nx.DiGraph()
        test_graph.add_edge(
            u="ENSDART00000149335.2:0-486",
            v="ENSDART00000149335.2:485-3379"
        )
        sealed_graph = correct_splice_graph(splice_graph, args)
        print(test_graph.nodes())
        print(test_graph.edges())
        print(sealed_graph.nodes())
        print(sealed_graph.edges())
        self.assertTrue(nx.is_isomorphic(
            sealed_graph,
            test_graph
        ))


if __name__ == '__main__':



    # write_gfa1(
    #     splice_graph=splice_graph,
    #     exons=exon_index,
    #     filename=args["output_gfa"]
    # )

    unittest.main()
    # Remove BF
    remove(temp_bloom, temp_gfa)
