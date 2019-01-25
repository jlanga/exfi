#!/usr/bin/env python3

"""
exfi.io.masking.py: submodule to soft and hard mask sequence strings
"""

import logging

import numpy as np

def cigar_to_int(cigar):
    """Convert a simple CIGAR string to overlap int

    >>> cigar_to_int('71N')
    -71
    >>> cigar_to_int('3M')
    3
    """
    if cigar[-1] == 'N':
        return -int(cigar[:-1])
    return int(cigar[:-1])


def soft_mask(sequence, left, right):
    """Lowercase the first left bases and last right bases of sequence

    >>> soft_mask('ACCGATCGATCGTAG', 2, 1)
    'acCGATCGATCGTAg'
    >>> soft_mask('ACCGATCGATCGTAG', 0, 2)
    'ACCGATCGATCGTag'
    >>> soft_mask('ACCGATCGATCGTAG', 2, 0)
    'acCGATCGATCGTAG'
    >>> soft_mask('ACCGATCGATCGTAG', 0, 0)
    'ACCGATCGATCGTAG'
    """
    if left == 0 and right == 0:
        return sequence
    if left == 0 and right > 0:
        return sequence[:-right] + sequence[-right:].lower()
    if left > 0 and right == 0:
        return sequence[:left].lower() + sequence[left:]
    return sequence[:left].lower() + sequence[left:-right] + sequence[-right:].lower()



def hard_mask(sequence, left, right):
    """Mask with N the first left bases and last right bases of sequence

    >>> hard_mask('ACCGATCGATCGTAG', 2, 1)
    'NNCGATCGATCGTAN'
    >>> hard_mask('ACCGATCGATCGTAG', 0, 2)
    'ACCGATCGATCGTNN'
    >>> hard_mask('ACCGATCGATCGTAG', 2, 0)
    'NNCGATCGATCGTAG'
    >>> hard_mask('ACCGATCGATCGTAG', 0, 0)
    'ACCGATCGATCGTAG'
    """
    if left == 0 and right == 0:
        return sequence
    if left == 0 and right > 0:
        return sequence[:-right] + 'N' * right
    if left > 0 and right == 0:
        return 'N' * left + sequence[left:]
    return 'N' * left + sequence[left:-right] + 'N' * right



def mask(node2sequence, edge2overlap, masking: str = "none"):
    """If any of the soft mask or hard mask are activated, mask

    :param dict exon_dict: Dict of the shape exon_id: sequence.
    :param dict overlap_dict: Dict of the shape (exon1, exon2): overlap between them.
    :param str masking: Type of masking to apply. Options: hard, soft, none
     (Default value = "None")    .
    """
    if masking == 'none':
        return node2sequence

    # Compose a dataframe of name, sequence, bases to trim to the left
    # and bases to trim to the right

    complete = node2sequence.merge(
        edge2overlap[['u', 'overlap']]\
            .rename(columns={'u': 'name', 'overlap': 'mask_right'}),
        on=['name'],
        how='outer'
    ).merge(
        edge2overlap[['v', 'overlap']]\
        .rename(columns={'v': 'name', 'overlap': 'mask_left'}),
        on=['name'],
        how='outer'
    )\
    .fillna(0)\
    .astype({'mask_right': np.int64, 'mask_left':np.int64})

    # Set to zero overlaps < 0
    complete['mask_right'] = complete.mask_right\
        .map(lambda x: x if x > 0 else 0)
    complete['mask_left'] = complete.mask_left\
        .map(lambda x: x if x > 0 else 0)

    if masking == "hard":
        logging.info("\tHard masking sequences")
        complete['sequence'] = complete.apply(
            lambda x: hard_mask(x.sequence, x.mask_left, x.mask_right),
            axis=1
        )
    elif masking == "soft":
        logging.info("\tSoft masking sequences")
        complete['sequence'] = complete.apply(
            lambda x: soft_mask(x.sequence, x.mask_left, x.mask_right),
            axis=1
        )

    node2sequence_masked = complete\
        [['name', 'sequence']]\
        .reset_index(drop=True)

    return node2sequence_masked
