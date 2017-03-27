#!/usr/bin/env python3

def trim_sequence(iterable_records, bases_to_the_left=0, bases_to_the_right=0):
    """(iterable, int, int) -> generator
    Given a record, trim the bases to the left and the ones to the right.
    """

    for record in iterable_records:
        if len(record) > bases_to_the_left + bases_to_the_right:
            n = len(record.seq)
            record.seq = record.seq[bases_to_the_left:n-bases_to_the_right]  
            yield record