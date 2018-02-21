#!/usr/bin/env python3

"""Module containing multiple functions for data wrantling:

- GFA1 <-> splice_graph
- GFA1 -> FASTA
- Soft/hard masking FASTA files
- bed3 to bed6
- Data frames, ...
"""

import logging

from typing import \
    Tuple

from exfi.classes import Coordinate


def _coordinate_tuple_to_str(coordinate: Coordinate) -> str:
    """Convert coordinates to str seqid:start-end

    :param str seqid: Sequence identifier
    :param int start: Start position
    :param int end: End position
    """
    seqid, start, end = coordinate
    return f"{seqid}:{start}-{end}"
