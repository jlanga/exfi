#!/usr/bin/env python3

from Bio.SeqIO.FastaIO import FastaWriter

def write_fasta(generator_seq, filename):
    """(filename, generator of Seqs) -> NoneType
    Write the fasta files from the generator into filename in fasta format 
    """
    with open(filename, "w") as handle_out:
        # SeqIO.write(
        #     format= "fasta",
        #     handle = handle_out,
        #     sequences = generator_seq
        # )
        writer = FastaWriter(handle_out, wrap= 0)
        writer.write_file(generator_seq)
