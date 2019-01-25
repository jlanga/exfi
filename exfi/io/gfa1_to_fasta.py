#!/usr/bin/env python3

"""exfi.io.gfa1_to_exons.py: submodule to read a gfa1, extract the exons and
store it in fasta format"""

from exfi.io.read_gfa import read_gfa1
from exfi.io.masking import mask, cigar_to_int

def gfa1_to_exons(fasta_out, gfa1_in, masking='none'):
    """Extract the exons in Fasta format"""
    with open(fasta_out, "w") as fasta:

        gfa1 = read_gfa1(gfa1_in)

        node2sequence = gfa1['segments']\
            .drop(columns='record_type')

        if node2sequence.shape[0] == 0:
            return

        edge2overlap = gfa1['links']\
            .drop(columns=['record_type', 'from_orient', 'to_orient'])\
            .rename(columns={
                'from': 'u', 'to': 'v', 'overlap': 'overlap_cigar'
            })
        edge2overlap['overlap'] = edge2overlap.overlap_cigar.map(cigar_to_int)

        node2sequence = mask(
            node2sequence=node2sequence,
            edge2overlap=edge2overlap,
            masking=masking
        )

        node2sequence["fasta"] = \
            ">" + node2sequence["name"] + "\n" + \
            node2sequence["sequence"]

        node2sequence.fasta.values.tofile(fasta, sep="\n", format="%s")
        fasta.write("\n")  # Final end line



def gfa1_to_gapped_transcripts(
        fasta_out, gfa1_in, gap_size=100, masking='none'):
    """Convert a GFA1 file to a gapped transcript file"""

    with open(fasta_out, "w") as fasta:

        separator = gap_size * 'N'

        gfa1 = read_gfa1(gfa1_in)

        # Segments -> node2sequence
        node2sequence = gfa1['segments']\
            .drop(columns=["record_type"])

        if node2sequence.shape[0] == 0:
            return

        edge2overlap = gfa1['links']\
            .drop(columns=['record_type', 'from_orient', 'to_orient'])\
            .rename(columns={
                'from': 'u', 'to': 'v', 'overlap': 'overlap_cigar'
            })
        edge2overlap["overlap"] = edge2overlap.overlap_cigar.map(cigar_to_int)

        path2nodes = gfa1['paths']\
            .drop(columns=['record_type', 'overlaps'])
        path2nodes.segment_names = path2nodes\
            .segment_names.str.replace('+', '')

        # Mask the sequences
        node2sequence = mask(
            node2sequence=node2sequence,
            edge2overlap=edge2overlap,
            masking=masking
        )

        node2sequence_dict = node2sequence\
            .set_index('name')\
            .to_dict()['sequence']

        # Compose the sequence
        path2nodes["gapped_sequence"] = path2nodes\
            .segment_names\
            .str.split(',')\
            .map(lambda x: separator.join([node2sequence_dict[y] for y in x]))

        # Create the fasta line
        path2nodes["fasta"] = \
            ">" + path2nodes.path_name + " " + path2nodes.segment_names + \
            "\n" + \
            path2nodes.gapped_sequence

        # Dump everything
        path2nodes.fasta.values.tofile(fasta, sep="\n", format="%s")
        fasta.write("\n")  # Final end line
