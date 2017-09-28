#!/usr/bin/env python3

from exfi.exons_to_gapped_transcript import \
    build_transcript_to_exon_dict
from Bio import SeqIO
import networkx as nx
import pandas as pd


def exons_to_df(exons):
    results = []
    for exon in exons.values():
        exon_id = exon.id
        transcript_coords = exon.description.split(" ")[1:]
        for transcript_coord in transcript_coords:
            transcript_id, coords = transcript_coord.split(":")
            start, end = coords.split("-")
            start = int(start)
            end = int(end)
            results.append(
                (transcript_id, start, end, exon_id, 0, "+")
            )
    return pd.DataFrame(
            data=results,
            columns=['transcript_id', 'start', 'end', 'exon_id', 'score', 'strand']
        )\
        .sort_values(
            by=['transcript_id', 'start','end']
        )
        
