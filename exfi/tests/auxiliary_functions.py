#!/usr/bin/env python3

"""
Auxiliary functions for testing
"""

def equal_seqrecord(seqrecord1, seqrecord2):
    """
    Check if seqrecord1 and seqrecord2 are identical in terms of:
    - id
    - sequence (as a Seq object)
    """
    if(seqrecord1.id == seqrecord2.id and seqrecord1.seq == seqrecord2.seq):
        return True
    else:
        return False


def equal_list_of_seqrecords(list_of_seqrecords1, list_of_seqrecords2):
    """
    Check if each element of list_of_seqrecords1 is exactly equal to each one of
    list_of_seqrecords2.
    """
    n1 = len(list_of_seqrecords1)
    n2 = len(list_of_seqrecords2)
    if n1 != n2:
        raise AssertionError(
            'Lengths differ:\n {len_1} != {len_2}'.format(
                len_1 = n1,
                len_2 = n2
            )
        )
    else:
        for i in range(n1):
            record1=list_of_seqrecords1[i]
            record2=list_of_seqrecords2[i]
            if not equal_seqrecord(record1, record2):
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
