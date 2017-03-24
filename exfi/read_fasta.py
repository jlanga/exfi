#!/usr/bin/env python3

from Bio import SeqIO

def read_fasta(fasta_filename):
    """(str) -> list of Seq
    Read a fasta file from fasta_filename and return each record via a 
    generator
    """
    with open(fasta_filename, "r") as handle_in:
        for record in SeqIO.parse(handle_in, "fasta"):
            yield record