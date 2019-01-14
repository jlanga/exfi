#!/usr/bin/env python3

"""
exfi.io.masking.py: submodule to soft and hard mask sequence strings
"""

import logging


def _process_overlap_cigar(cigar_string: str) -> list:
    """Process a simple CIGAR string.

    :param cigar_string: Process a simple CIGAR string of the shape number-letter: 90G, 10M, ...
    """
    return [cigar_string[-1], int(cigar_string[:-1])]


def _soft_mask_right(string: str, n_bases: int) -> str:
    """Soft mask the rightmost n bases.

    :param str string: string of nucleotides to be soft masked.
    :param int n_bases: number of letters at the right end (3') to be soft masked.
    """
    return string[:-n_bases] + string[-n_bases:].lower()


def _soft_mask_left(string: str, n_bases: int) -> str:
    """Soft mask the leftmost n bases.

    :param str string: string of nucleotides to be soft masked
    :param int n_bases: number of letters at the left end (5') to be soft masked.
    """
    return string[:n_bases].lower() + string[n_bases:]


def _soft_mask(exon_dict, overlap_dict):
    """Soft mask all overlaps in the exon_dict.

    :param dict exon_dict: dict of exon_id: sequence
    :param dict overlap_dict: dict of (node1, node2): overlap
    """
    exon_dict = exon_dict.copy()
    for (start, end), overlap in overlap_dict.items():
        if overlap > 0:
            exon_dict[start] = _soft_mask_right(exon_dict[start], overlap)
            exon_dict[end] = _soft_mask_left(exon_dict[end], overlap)
    return exon_dict


def _hard_mask_right(string: str, n_bases: int) -> str:
    """Hard mask the rightmost n_bases bases

    :param string: Nucleotide sequence to hard mask.
    :param n_bases: Number of bases to hard mask at the right (3') end.
    """
    return string[:-n_bases] + "N" * n_bases


def _hard_mask_left(string: str, n_bases: int):
    """Hard mask the leftmost n_bases bases

    :param str string: Nucleotide sequence to hard mask.
    :param int n_bases: Number of bases to hard mask at the left (5') end.
    """
    return "N" * n_bases + string[n_bases:]


def _hard_mask(exon_dict, overlap_dict):
    """Hard mask all overlaps in the exon_dict.

    :param dict exon_dict: Dict of the shape exon_id: sequence.
    :param dict overlap_dict: Dict of the shape (exon1, exon2): overlap between them.
    """
    exon_dict = exon_dict.copy()
    for (start, end), overlap in overlap_dict.items():
        if overlap > 0:
            exon_dict[start] = _hard_mask_right(exon_dict[start], overlap)
            exon_dict[end] = _hard_mask_left(exon_dict[end], overlap)
    return exon_dict


def _mask(exon_dict, overlap_dict, masking: str = "none"):
    """If any of the soft mask or hard mask are activated, mask

    :param dict exon_dict: Dict of the shape exon_id: sequence.
    :param dict overlap_dict: Dict of the shape (exon1, exon2): overlap between them.
    :param str masking: Type of masking to apply. Options: hard, soft, none (Default value = "None")
    .
    """
    if masking == "hard":
        logging.info("\tHard masking sequences")
        exon_dict = _hard_mask(exon_dict, overlap_dict)
    elif masking == "soft":
        logging.info("\tSoft masking sequences")
        exon_dict = _soft_mask(exon_dict, overlap_dict)
    return exon_dict
