#!/usr/bin/env python3

"""exfi.io.fasta_to_dict.py: submodule to convert a fasta file into a dict"""

from Bio.SeqIO.FastaIO import \
    SimpleFastaParser

def fasta_to_dict(filename: str) -> dict:
    """Fast Fasta to dict via SimpleFastaParser

    :param filename: str: Path to the fasta file

    """
    with open(filename, "r") as handle:
        return dict(SimpleFastaParser(handle))
