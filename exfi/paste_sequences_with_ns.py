#!/usr/bin/env python3

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def paste_sequences_with_ns(iterable_records, number_of_ns=100):
    """(iterable, int) -> generator
    Join sequences from iterable_records with Ns
    """

    current_record = None

    for record in iterable_records:

        record.id = record.id.rsplit("_e")[0]

        if current_record and record.id == current_record.id:
            current_record.seq = current_record.seq + "N" * number_of_ns + str(record.seq)
        else:
            yield current_record
            current_record = record