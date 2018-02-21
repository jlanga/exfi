#!/usr/bin/env python3

"""
tests.custom_assertions.py: custom assertions for unit tests:
- assertEqualListOfSeqrecords: check if a list of seqrecords have:
    - the same length
    - the same id
    - the same sequence
- assertEqualSpliceGraphs: check if two splice graphs:
    - are isomorphic with nx.is_isomorphic
    - each node have the same coordinates
    - each edge have the same overlap
"""

from typing import List, Dict

import networkx as nx

import pandas as pd

from Bio.SeqRecord import SeqRecord

from exfi.classes import SpliceGraph, SpliceGraphDict

def check_same_keys(dict1: dict, dict2: dict) -> None:
    """Check if two dicts have the exact same keys"""
    if set(dict1.keys()) != set(dict2.keys()):
        raise KeyError("Keys differ: {keys1} {keys2}".format(
            keys1=dict1.keys(), keys2=dict2.keys()
        ))



def check_same_values(dict1: dict, dict2: dict) -> None:
    """Check if two dicts have the same values"""
    for key, value1 in dict1.items():  # Check same values
        value2 = dict2[key]
        if value1 != value2:
            raise ValueError("{key1}: {value1} != {key2} : {value2}".format(
                key1=key, value1=value1, key2=key, value2=value2
            ))



def check_same_dict(dict1: dict, dict2: dict) -> None:
    """Check if two dicts contain the exact same values"""
    check_same_keys(dict1, dict2)
    check_same_values(dict1, dict2)



def check_equal_node2coord(sg1: SpliceGraph, sg2: SpliceGraph) -> None:
    """Check if two splice graphs have the same node2coord dicts"""
    node2coord1 = nx.get_node_attributes(G=sg1, name="coordinates")
    node2coord2 = nx.get_node_attributes(G=sg2, name="coordinates")
    check_same_dict(node2coord1, node2coord2)



def check_equal_edge2overlap(sg1: SpliceGraph, sg2: SpliceGraph) -> None:
    """Check if two splice graphs have the same node2coord dicts"""
    edge2overlap1 = nx.get_edge_attributes(G=sg1, name="overlaps")
    edge2overlap2 = nx.get_edge_attributes(G=sg2, name="overlaps")
    check_same_dict(edge2overlap1, edge2overlap2)



def check_equal_df_dict_values(dict1: dict, dict2: dict) -> None:
    """Check if two data frames are equal
    Solution: https://stackoverflow.com/a/33223893
    """
    from numpy import array_equal
    for key, df1 in dict1.items():
        df2 = dict2[key]
        if not array_equal(df1, df2):
            raise ValueError("df1 != df2:\n{df1}\n{df2}".format(df1=df1, df2=df2))



def check_equal_splice_graphs(sg1: SpliceGraph, sg2: SpliceGraph) -> None:
    """Check if two splice graphs are:
        - isomorphic
        - node2coord are equal
        - edge2overlaps are equal
    """
    if not nx.is_isomorphic(sg1, sg2):
        AssertionError("splicegraph are not isomorphic")
    check_equal_node2coord(sg1, sg2)
    check_equal_edge2overlap(sg1, sg2)



def check_equal_dict_of_sg(dict1: dict, dict2: dict) -> None:
    """Check if each key, element are equal splice graphs"""
    check_same_keys(dict1, dict2)
    for key, sg1 in dict1.items():
        sg2 = dict2[key]
        check_equal_splice_graphs(sg1, sg2)



def check_equal_length(iter1: List, iter2: List) -> None:
    """Check if two iterables have the same length"""
    length_1 = len(iter1)
    length_2 = len(iter2)
    if length_1 != length_2:
        raise AssertionError('Lengths differ: {len_1} != {len_2}'.format(
            len_1=length_1, len_2=length_2
        ))



def check_equal_seqrecrods(seqrecord1: SeqRecord, seqrecord2: SeqRecord) -> None:
    """Check if id and seq are equal"""
    if seqrecord1.id != seqrecord2.id or seqrecord1.seq != seqrecord2.seq:
        raise AssertionError(
            'Records differ: {id1}: {seq1} {id2}: {seq2}'.format(
                id1=seqrecord1.id, seq1=seqrecord1.seq, id2=seqrecord2.id, seq2=seqrecord2.seq
            )
        )



def check_equal_list_seqrecords(iter1: List[SeqRecord], iter2: List[SeqRecord]) -> None:
    """Check if a list of SeqRecords are equal"""
    for i, _ in enumerate(iter1):
        check_equal_seqrecrods(iter1[i], iter2[i])



class CustomAssertions:
    """
    Custom assertions not covered in unittest:
    - assertEqualListOfSeqrecords
    """

    @classmethod
    def assertEqualDict(self, dict1: dict, dict2: dict) -> None:
        """Check if two dicts are equal (values are compared with ==)"""
        # pylint: disable=invalid-name, bad-classmethod-argument
        check_same_dict(dict1, dict2)



    @classmethod
    def assertEqualListOfSeqrecords(
            self, records1: List[SeqRecord], records2: List[SeqRecord]) -> None:
        """
        Check if each element of list_of_seqrecords1 is exactly equal to each one of
        list_of_seqrecords2.
        """
        # pylint: disable=invalid-name, bad-classmethod-argument
        check_equal_length(records1, records2)
        check_equal_list_seqrecords(records1, records2)



    @classmethod
    def assertEqualSpliceGraphs(self, sg1: SpliceGraph, sg2: SpliceGraph) -> None:
        """Check if two splice graph are equal:"""
        # pylint: disable=invalid-name,bad-classmethod-argument
        check_equal_splice_graphs(sg1, sg2)



    @classmethod
    def assertEqualDictOfDF(
            self, dict1: Dict[str, pd.DataFrame], dict2: Dict[str, pd.DataFrame]) -> None:
        """Check if two dicts of pd.DataFrame are equal"""
        # pylint: disable=invalid-name,bad-classmethod-argument
        check_same_keys(dict1, dict2)
        check_equal_df_dict_values(dict1, dict2)


    @classmethod
    def assertEqualDictOfSpliceGraphs(self, dict1: SpliceGraphDict, dict2: SpliceGraphDict) -> None:
        """Check if two dicts of nx.DiGraph and some data attached to nodes and edges are equal"""
        # pylint: disable=invalid-name, bad-classmethod-argument
        check_equal_dict_of_sg(dict1, dict2)
