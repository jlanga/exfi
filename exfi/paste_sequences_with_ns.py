#!/usr/bin/env python3

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def paste_sequences_with_ns(iterable_records, number_of_ns=100):
    """(iterable, int) -> generator
    Join sequences from iterable_records with Ns
    """

    # store and process the transcripts into a dict record.id - seq
    transcript_to_seq = {}

    for record in iterable_records:

        record.id = record.id.rsplit("_e")[0]
         
        if record.id not in transcript_to_seq:
            transcript_to_seq[record.id] = str(record.seq)
        else:
            transcript_to_seq[record.id] += "N" * number_of_ns + str(record.seq)
    
    for record_name, sequence in transcript_to_seq.items():
        yield SeqRecord(
            id = record_name,
            seq = Seq(sequence)
        )