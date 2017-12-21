from unittest import TestCase

from exfi.correct_splicegraph import \
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

from exfi.build_splicegraph import \
    build_splicegraph

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
    "input_fasta": "tests/files/correct_splicegraph/transcript.fa",
    "max_fp_bases": 5,
    "max_overlap": 10,
    "output_gfa": temp_gfa[1],
    "threads": 4,
    "max_gap_size": 10,
    "reads": ["tests/files/correct_splicegraph/reads.fa"]
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
positive_exons_fasta = _get_fasta(transcriptome_index, positive_exons_bed)

# Reduce and convert to dict
# positive_exons_fasta = reduce_exons(positive_exons_fasta)
exon_index = {exon.id: exon for exon in positive_exons_fasta}
# reduced_exons = reduce_exons(positive_exons_fasta)
# exon_index = {exon.id: exon for exon in reduced_exons}

# Build splice graph
splice_graph = build_splicegraph(exon_index)




class TestCoordinatesToVariables(TestCase):
    """ _coordinates_to_variables(coordinates):
    (string) -> (string, int, int)
    """
    def test_empty_string(self):
        """Split an empty string"""
        with self.assertRaises(ValueError):
            _coordinates_to_variables("")

    def test_incorrect_string(self):
        """Split a string not in the form id:start-end"""
        with self.assertRaises(IndexError):
            _coordinates_to_variables("ENSDART00000161035:1")

    def test_correct_string(self):
        """Split a standard string"""
        self.assertEqual(
            _coordinates_to_variables("ENSDART00000161035:1-15"),
            ("ENSDART00000161035", 1, 15)
        )

    def test_messy_string(self):
        """Split a string with multiple _ and :"""
        self.assertEqual(
            _coordinates_to_variables("TRINITY_g14_c15_i5:1:2-15-12:1-15"),
            ("TRINITY_g14_c15_i5:1:2-15-12", 1, 15)
        )



class TestPrepareSealer(TestCase, CustomAssertions):
    """_prepare_sealer(splicegraph, args)
    (nx.DiGraph, dict_of_parameters) -> str
    """
    def test_file_creation(self):
        """Create the fasta file for sealer"""
        sealer_input_fn = _prepare_sealer(splice_graph, args)
        self.assertEqualListOfSeqrecords(
            list(SeqIO.parse(
                sealer_input_fn,
                format="fasta"
            )),
            list(SeqIO.parse(
                "tests/files/correct_splicegraph/to_seal.fa",
                format="fasta"
            ))
        )
        remove(sealer_input_fn)



class TestRunSealer(TestCase, CustomAssertions):
    """_run_sealer(sealer_input_fn, args):
    (str, dict) -> str
    """
    def test_run(self):
        """Run sealer"""
        sealer_in_fn = _prepare_sealer(splice_graph, args)
        sealer_out_fn = _run_sealer(
            sealer_input_fn=sealer_in_fn,
            args=args
        )
        self.assertEqualListOfSeqrecords(
            list(SeqIO.parse(
                "tests/files/correct_splicegraph/sealed.fa",
                format="fasta"
            )),
            list(SeqIO.parse(
                sealer_out_fn, "fasta"
            ))
        )
        remove(sealer_input_fn)
        remove(sealer_output_fn)



class TestCollectSealerResults(TestCase):
    """_collect_sealer_results(handle):
    (str) -> dict
    """
    def test_collect_empty(self):
        pass

    def test_collect_somedata(self):
        pass



class TestSculptGraph(TestCase):
    """_sculpt_graph(splice_graph, edge2fill):
    (nx.DiGraph, dict) -> nx.DiGraph
    """
    def test_sculpt_empty_data(self):
        pass
    def test_sculpt_wrong_data(self):
        pass
    def test_scuplt_real_data(self):
        pass



class TestCorrectSpliceGraph(TestCase):
    """ correct_splice_graph(splice_graph, args):
    (nx.DiGraph, int) -> nx.DiGraph
    """
    def test_correct_splice_graph(self):
        pass


if __name__ == '__main__':



    # write_gfa1(
    #     splice_graph=splice_graph,
    #     exons=exon_index,
    #     filename=args["output_gfa"]
    # )

    unittest.main()
    # Remove BF
    remove(temp_bloom, temp_gfa)
