#!/usr/bin/env python3

"""Module containing multiple functions for data wrantling:

- GFA1 <-> splice_graph
- GFA1 -> FASTA
- Soft/hard masking FASTA files
- bed3 to bed6
- Data frames, ...
"""

import logging


def _coordinate_str_to_tuple(coordinates: str) -> tuple:
    """Convert a string in ":,-" format into a tuple

    :param str coordinates: string of the form "seq:start-end"
    """
    transcript, start_end = coordinates.rsplit(":", 1)
    start_end = start_end.rsplit("-", 1)
    start = int(start_end[0])
    end = int(start_end[1])
    return transcript, start, end


def _coordinate_tuple_to_str(seqid: str, start: int, end: int) -> str:
    """Convert coordinates to str seqid:start-end

    :param str seqid: Sequence identifier
    :param int start: Start position
    :param int end: End position
    """
    return f"{seqid}:{start}-{end}"
