#!/usr/bin/env python3

"""exfi.io.fasta_to_dict.py: submodule to convert a fasta file into a dict"""

import logging

from Bio.SeqIO.FastaIO import \
    SimpleFastaParser

# from exfi.classes import FastaDict

def fasta_to_dict(filename):
    """Fast Fasta to dict via SimpleFastaParser

    :param filename: str: Path to the fasta file

    """
    logging.info('Dumping fasta to dict')
    with open(filename, "r") as handle:
        return {
            identifier.split()[0]: sequence
            for identifier, sequence in SimpleFastaParser(handle)
        }
    logging.info('Done')
