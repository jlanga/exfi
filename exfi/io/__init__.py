#!/usr/bin/env python3

"""
Module containing multiple functions for data wrantling:
GFA1 <-> splice_graph
bed3 to bed6
Data frames, ...
"""


from Bio import SeqIO


def _coordinate_str_to_tuple(coordinates):
    """(string) -> (string, int, int)

    Convert a genomic coordinate of the form "TR_ID:start-end" into the tuple
    ("TR_ID", start, end).
    """
    transcript, start_end = coordinates.rsplit(":", 1)
    start_end = start_end.rsplit("-", 1)
    start = int(start_end[0])
    end = int(start_end[1])
    return transcript, start, end


def _coordinate_tuple_to_str(transcript_id, start, end):
    """(str, int, int) -> str

    Convert coordinates to str transcript_id:start-end
    """
    return "{transcript_id}:{start}-{end}".format(
        transcript_id=transcript_id,
        start=start,
        end=end
    )

def _clean_seqrecord(seqrecord):
    """Delete the identifier from the description"""
    seqrecord.description = " ".join(seqrecord.description.split(" ")[1:])
    return seqrecord



def _clean_index(index):
    """Clean all elements from an indexed fasta"""
    index_clean = {}
    for key, value in index.items():
        index_clean[key] = _clean_seqrecord(value)
    return index_clean



def index_fasta(filename):
    """Create a fasta dict, with clean descriptions, key=id, value=seqrecord"""
    index = SeqIO.index(
        filename=filename,
        format="fasta"
    )
    return _clean_index(index)
