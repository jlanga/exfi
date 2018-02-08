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



class CustomAssertions:
    """
    Custom assertions not covered in unittest:
    - assertEqualListOfSeqrecords
    """
    @classmethod
    def assertEqualListOfSeqrecords(self, records1, records2):
        """
        Check if each element of list_of_seqrecords1 is exactly equal to each one of
        list_of_seqrecords2.
        """
        # pylint: disable=invalid-name, bad-classmethod-argument
        length_1 = len(records1)
        length_2 = len(records2)
        if length_1 != length_2:
            raise AssertionError(
                'Lengths differ:\n {len_1} != {len_2}'.format(
                    len_1=length_1,
                    len_2=length_2
                )
            )
        else:
            for i in range(length_1):
                record1 = records1[i]
                record2 = records2[i]
                if record1.id != record2.id or record1.seq != record2.seq:
                    raise AssertionError(
                        'Records at position {i} differ:\n{id1} : {seq1}\n{id2} : {seq2}'.format(
                            i=i,
                            id1=record1.id,
                            seq1=record1.seq,
                            id2=record2.id,
                            seq2=record2.seq
                        )
                    )
            return True

    @classmethod
    def assertEqualSpliceGraphs(self, sg1, sg2):
        """
        Check if two splice graph are equal:
        - are isomorphic
        - same coordinates
        - same overlaps
        """
        # pylint: disable=invalid-name,bad-classmethod-argument

        import networkx as nx

        coordinates1 = nx.get_node_attributes(G=sg1, name="coordinates")
        coordinates2 = nx.get_node_attributes(G=sg1, name="coordinates")
        overlaps1 = nx.get_edge_attributes(G=sg1, name="overlaps")
        overlaps2 = nx.get_edge_attributes(G=sg1, name="overlaps")
        same_coordinates = coordinates1 == coordinates2
        same_overlaps = overlaps1 == overlaps2
        are_isomorphic = nx.is_isomorphic(sg1, sg2)
        if are_isomorphic and same_coordinates and same_overlaps:
            return True
        return False

    @classmethod
    def assertEqualDictOfDF(self, dict1, dict2):
        """Check if two dicts of pd.DataFrame are equal"""
        # pylint: disable=invalid-name,bad-classmethod-argument

        for key, dataframe in dict1.items():

            if key not in dict2:
                raise KeyError("Key %s not in dict2" % key)
            else:
                if dataframe.equals(dict2[key]):
                    return True
                return False

    @classmethod
    def assertEqualDictOfSpliceGraphs(self, dict1, dict2):
        """Check if two dicts of nx.DiGraph and some data attached to nodes and edges are equal"""
        # pylint: disable=invalid-name, bad-classmethod-argument


        if len(dict1) != len(dict2):
            raise AssertionError("lengths differ!")

        if not dict1:
            return True

        for key, sg1 in dict1.items():
            if key not in dict2:
                raise KeyError("Key %s not in dict2" % key)
            else:
                sg2 = dict2[key]

                import networkx as nx

                coordinates1 = nx.get_node_attributes(G=sg1, name="coordinates")
                coordinates2 = nx.get_node_attributes(G=sg1, name="coordinates")
                overlaps1 = nx.get_edge_attributes(G=sg1, name="overlaps")
                overlaps2 = nx.get_edge_attributes(G=sg1, name="overlaps")
                same_coordinates = coordinates1 == coordinates2
                same_overlaps = overlaps1 == overlaps2
                are_isomorphic = nx.is_isomorphic(sg1, sg2)
                if are_isomorphic and same_coordinates and same_overlaps:
                    return True
                return False
