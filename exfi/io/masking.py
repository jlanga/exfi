#!/usr/bin/env python3

"""
exfi.io.masking.py: submodule to soft and hard mask sequence strings
"""

import logging

import pandas as pd

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
    return  sequence[:left].lower() + sequence[left:-right] + sequence[-right:].lower()



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

    edge2overlap['tmp_overlap'] = edge2overlap.overlap.map(
        lambda x: x if x > 0 else 0
    )

    tmp = pd.merge(
        node2sequence,
        edge2overlap[['u', 'tmp_overlap']].rename(
            columns={'u': 'name', 'tmp_overlap': 'mask_right'}
        ),
        on=['name']
    )

    complete = complete = pd.merge(
        tmp,
        edge2overlap[['v', 'tmp_overlap']].rename(
            columns={'v': 'name', 'tmp_overlap': 'mask_left'}
        ),
        on=['name']
        )

    complete['tmp'] = tuple(zip(
        complete.sequence, complete.mask_left, complete.mask_right
    ))



    if masking == "hard":
        logging.info("\tHard masking sequences")
        complete['sequence'] = complete.tmp.map(lambda x: hard_mask(*x))
    elif masking == "soft":
        logging.info("\tSoft masking sequences")
        complete['sequence'] = complete.tmp.map(lambda x: soft_mask(*x))

    exon_dict = complete[['name', 'sequence']].reset_index(drop=True)
    return exon_dict
