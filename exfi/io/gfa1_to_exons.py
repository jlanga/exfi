#!/usr/bin/env python3

"""
exfi.io.gfa1_to_exons.py: submodule to convert a GFA1 file into a fasta containing only the exons
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from exfi.io import _coordinate_tuple_to_str
from exfi.io.read_gfa1 import read_gfa1
from exfi.io.masking import _mask



def gfa1_to_exons(gfa_in_fn, fasta_out_fn, soft_mask_overlaps=False, hard_mask_overlaps=False):
    """(str, str, bool, bool) -> None

    Write the exons in FASTA format present in a GFA1 file
    """
    gfa1 = read_gfa1(gfa_in_fn)

    exon2sequence = gfa1["segments"]
    exon2coordinates = gfa1["containments"]
    link2overlap = gfa1["links"]

    # Mask if necessary
    exon2sequence = _mask(exon2sequence, link2overlap, soft_mask_overlaps, hard_mask_overlaps)

    # Add coordinate information to description
    # Compose SeqRecord of each exon
    sequences = []
    for exon_id, exon_sequence in exon2sequence.items():
        # Compose coordinates
        exon_coordinates = exon2coordinates[exon_id]
        description = " ".join(
            _coordinate_tuple_to_str(*coordinate)
            for coordinate in exon_coordinates
        )
        sequences.append(SeqRecord(
            id=exon_id,
            description=description,
            seq=Seq(exon_sequence)
        ))

    # Write to fasta
    SeqIO.write(format="fasta", handle=fasta_out_fn, sequences=sequences)
