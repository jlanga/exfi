#!/usr/bin/env python3

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def paste_sequences_with_ns(iterable_records, number_of_ns=100):
    """(iterable, int) -> generator
    Join sequences from iterable_records with Ns
    """
    ns = "N" * number_of_ns

    current_id = None
    current_seq = []

    for record in iterable_records:

        record.id = record.id.rsplit("_e")[0]

        if current_id and current_id != record.id:
            yield SeqRecord(id = current_id, seq = Seq(ns.join(current_seq)), description = "")
            current_id = None
            current_seq = []
        
        current_id = record.id
        current_seq.append(str(record.seq))

    if current_id:
        yield SeqRecord(id = current_id, seq = Seq(ns.join(current_seq)), description = "")


