#!/usr/bin/env python3


# Imports
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def reduce_exons(iterable_of_seqrecords):
    """Reduce the exons by sequence identity

    Build 2 dicts in which keys are the sequences. In one of them store
    {seq : exon_identifier} and in the other {seq : description} where
    description is the exon coordinates.
    Once collected all records, return them.
    """
    # Initialize
    seq_to_id = dict()
    seq_to_description = dict()

    # Go over each record
    for record in iterable_of_seqrecords:
        seq = str(record.seq)
        if seq in seq_to_id.keys():  # Just append
            seq_to_description[seq] += " " + record.description
        else:  # Enter values in both dicts
            seq_to_description[seq] = record.description
            seq_to_id[seq] = 'EXON{:011d}'.format(len(seq_to_id) + 1)

    # Collect data from both dicts into a SeqRecord and return
    for seq in seq_to_id.keys():
        record = SeqRecord(
            id=seq_to_id[seq],
            description=seq_to_description[seq],
            seq=Seq(seq)
        )
        yield record
