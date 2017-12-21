from unittest import TestCase
from exfi.correct_splicegraph import \
    _coordinates_to_variables, \
    _prepare_sealer, \
    _run_sealer, \
    _collect_sealer_results , \
    _sculpt_graph, \
    correct_splice_graph

import tempfile


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



class TestPrepareSealer(TestCase):
    """_prepare_sealer(splicegraph, args)
    (nx.DiGraph, dict_of_parameters) -> str
    """
    def test_file_creation(self):
        pass



class TestRunSealer(TestCase):
    """_run_sealer(sealer_input_fn, args):
    (str, dict) -> str
    """
    def test_run(self):
        pass



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
