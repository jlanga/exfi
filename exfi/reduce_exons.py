#!/usr/bin/env python3


# Imports
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def reduce_exons(iterable_of_seqrecords):
    """Reduce the exons by sequence identity

    Build a dict whose keys are the sequence in str format and the values are
    lists in which the first element is the exon_id and the remaining are
    the coordinates in loci:start-end format
    """
    # Initialize
    seq_to_data = dict()

    # Go over each record
    for record in iterable_of_seqrecords:
        seq = str(record.seq)
        if seq in seq_to_data:  # Just append
            seq_to_data[seq].append(record.description)
        else:  # Enter values in both dicts
            seq_to_data[seq] = [
                'EXON{:011d}'.format(len(seq_to_data) + 1),
                record.description
            ]

    # Collect data from both dicts into a SeqRecord and return
    for seq in seq_to_data.keys():
        identifier, *description = seq_to_data[seq]
        record = SeqRecord(
            id=identifier,
            description=" ".join(description),
            seq=Seq(seq)
        )
        yield record
