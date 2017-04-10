#!/usr/bin/env python3

from Bio import SeqIO

def read_fasta_from_bedtools_getfasta(iterable):

    current_transcript_id = ""
    exon_number = 0

    for record in iterable:
        # Get just the transcript id, not the coordinates
        transcript_id = record.id 
        
        # If we are looking at a different transcript
        if current_transcript_id != transcript_id:
            current_transcript_id = transcript_id
            record.id = record.description + ";1"
            record.description = ""
            yield record
            exon_number = 1
        else:
            exon_number += 1
            record.id = record.description + ";{}".format(exon_number)
            record.description = ""
            yield record
