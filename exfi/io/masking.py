#!/usr/bin/env python3

"""
exfi.io.masking.py: submodule to soft and hard mask sequence strings
"""

def _process_overlap_cigar(cigar_string):
    """Process a simple CIGAR string (number, letter)"""
    return [cigar_string[-1], int(cigar_string[:-1])]



def _soft_mask_right(string, n_bases):
    """Soft mask the rightmost n bases"""
    return string[:-n_bases] + string[-n_bases:].lower()



def _soft_mask_left(string, n_bases):
    """Soft mask the leftmost n bases"""
    return string[:n_bases].lower() + string[n_bases:]



def _soft_mask(exon_dict, overlap_dict):
    """Soft mask all overlaps in the exon_dict"""
    exon_dict = exon_dict.copy()
    for (start, end), overlap in overlap_dict.items():
        if overlap > 0:
            exon_dict[start] = _soft_mask_right(exon_dict[start], overlap)
            exon_dict[end] = _soft_mask_left(exon_dict[end], overlap)
    return exon_dict



def _hard_mask_right(string, n_bases):
    """Hard mask the rightmost n bases"""
    return string[:-n_bases] + "N" * n_bases



def _hard_mask_left(string, n_bases):
    """Hard mask the leftmost n bases"""
    return "N" * n_bases + string[n_bases:]



def _hard_mask(exon_dict, overlap_dict):
    """Hard mask all overlaps in the exon_dict"""
    exon_dict = exon_dict.copy()
    for (start, end), overlap in overlap_dict.items():
        if overlap > 0:
            exon_dict[start] = _hard_mask_right(exon_dict[start], overlap)
            exon_dict[end] = _hard_mask_left(exon_dict[end], overlap)
    return exon_dict



def _mask(exon_dict, overlap_dict, soft_mask_overlaps=False, hard_mask_overlaps=False):
    """If any of the soft mask or hard mask are activated, mask"""
    if soft_mask_overlaps and hard_mask_overlaps:
        raise Exception("I can't soft mask and hard mask at the same time, dude!")
    if soft_mask_overlaps:
        exon_dict = _soft_mask(exon_dict, overlap_dict)
    if hard_mask_overlaps:
        exon_dict = _hard_mask(exon_dict, overlap_dict)
    return exon_dict
