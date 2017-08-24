#!/usr/bin/env python3

"""
Auxiliary functions and classes for testing
"""

from Bio import SeqIO

class CustomAssertions:

    def assertEqualListOfSeqrecords(self, records1, records2):
        """
        Check if each element of list_of_seqrecords1 is exactly equal to each one of
        list_of_seqrecords2.
        """
        n1 = len(records1)
        n2 = len(records2)
        if n1 != n2:
            raise AssertionError(
                'Lengths differ:\n {len_1} != {len_2}'.format(
                    len_1 = n1,
                    len_2 = n2
                )
            )
        else:
            for i in range(n1):
                record1=records1[i]
                record2=records2[i]
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
