#!/usr/bin/env python3

"""exfi.io.gff3_to_bed.py: exfi submodule to convert a gff3 to bed3 where
coordinates are with respect to the transcriptome"""

import sys

import pandas as pd
import numpy as np

def gff3_to_bed3(gff3_in, mode="ensembl"):
    """Read a GFF3 file and convert it to BED3, where coordinates are with
    respect to the transcriptome

    Modes available:
        - "ensembl": for files downloaded from Ensembl,
        - "gmap": for GFF3 files generated from GMAP,
        - "ncbi": for GFF3 files downloaded from NCBI Genomes
    """

    gff3_columns = [
        "seqid", "source", "type", "start", "end", "score", "strand", "phase",
        "attributes"
    ]

    bed3_columns = ["chrom", "chromStart", "chromEnd"]
    bed3_dtypes = {
        "chrom": np.str, "chromStart": np.int64, "chromEnd": np.int64
    }

    raw = pd.read_csv(
        sep='\t',
        na_values=".",
        usecols=["type", "start", "end", "strand", "attributes"],
        filepath_or_buffer=gff3_in,
        comment="#",
        header=None,
        names=gff3_columns,
        low_memory=False  # Convert types at the end. Seqid is char, not int
    )

    if raw.shape[0] == 0:
        exons = pd.DataFrame(columns=bed3_columns)
        exons = exons.astype(bed3_dtypes)
        return exons

    if mode == "gmap":
        exons = raw[raw['type'] == 'cDNA_match'].drop(columns='type')
        exons['transcript_id'] = exons['attributes']\
            .str.split(";").str[1]\
            .str.extract(r'Name=([\w\d.-_]+)')
    elif mode == "ensembl":
        exons = raw[raw['type'] == 'exon'].drop(columns='type')
        exons["transcript_id"] = exons["attributes"]\
            .str.split(";", 1, ).str[0]\
            .str.extract(r'Parent=transcript:([\w\d.-_]+)')
    else:
        sys.exit("Unknown mode")


    if exons.shape[0] == 0:
        return exons

    exons = exons[['transcript_id', 'strand', 'start', 'end']]

    positive = (
        exons
        [exons['strand'] == '+']
        .drop(columns='strand')
        .sort_values(by=['transcript_id', 'start', 'end'])
    )

    negative = (
        exons
        [exons['strand'] == '-']
        .drop(columns='strand')
        .sort_values(
            by=['transcript_id', 'start', 'end'],
            ascending=[True, False, False]
        )
    )

    merged = pd.concat([positive, negative])

    merged['length'] = merged['end'] - merged['start'] + 1
    merged['transcript_end'] = (
        merged
        .groupby('transcript_id')
        ['transcript_id', 'length']
        .cumsum()
    )

    merged['transcript_start'] = merged['transcript_end'] - merged['length']

    merged = merged[['transcript_id', 'transcript_start', 'transcript_end']]

    merged = merged.rename(columns={
        'transcript_id': 'chrom',
        'transcript_start': 'chromStart',
        'transcript_end': 'chromEnd'
    })

    merged = merged.astype(bed3_dtypes)

    merged = merged.reset_index(drop=True)

    return merged
