#!/usr/bin/env python3

"""exfi.io.fasta_to_dict.py: submodule to convert a fasta file into a dict"""

from Bio.SeqIO.FastaIO import \
    SimpleFastaParser

from exfi.classes import FastaDict

def fasta_to_dict(filename: str) -> FastaDict:
    """Fast Fasta to dict via SimpleFastaParser

    :param filename: str: Path to the fasta file

    """
    with open(filename, "r") as handle:
        return FastaDict(
            (identifier.split()[0], sequence)
            for identifier, sequence in SimpleFastaParser(handle)
        )
