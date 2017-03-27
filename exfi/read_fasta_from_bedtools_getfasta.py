#!/usr/bin/env python3

from Bio import SeqIO

def read_fasta_from_bedtools_getfasta(iterable):

    current_transcript_id = ""
    exon_number = 0

    for record in iterable:
        record.id = record.id.rsplit(":")[0]
        if current_transcript_id != record.id:
            current_transcript_id = record.id
            record.id = record.id + "_e1"
            yield record
            exon_number = 1
        else:
            exon_number += 1
            record.id = record.id + "_e{}".format(exon_number)
            yield record
