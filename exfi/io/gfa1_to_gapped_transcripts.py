#!/usr/bin/env python3

"""
exfi.io.gfa1_to_gapped_transcript.py: submodule to convert a GFA1 file to a fasta file where the
spaces between predicted exons are filled with a string of Ns.
"""

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from exfi.io.masking import _mask
from exfi.io.read_gfa1 import read_gfa1


def _compose_paths(exon_dict, path_dict, number_of_ns):
    """(dict, dict, int) -> iterable of SeqRecord

    Compose and return each gapped transcript.
    """
    chunk_of_ns = "N" * number_of_ns
    for transcript_id, exon_list in sorted(path_dict.items()):
        exon_seqs = [str(exon_dict[exon_id]) for exon_id in exon_list]
        yield SeqRecord(
            id=transcript_id,
            description=",".join(exon_list),
            seq=Seq(chunk_of_ns.join(exon_seqs))
        )


def gfa1_to_gapped_transcripts(gfa_in, fasta_out, number_of_ns=100, masking="none"):
    """
    Write gapped transcripts as fasta from GFA1 file
    """

    # Process
    gfa1 = read_gfa1(gfa_in)
    exon_dict = gfa1["segments"]
    overlap_dict = gfa1["links"]
    path_dict = gfa1["paths"]


    # Mask if necessary
    exon_dict = _mask(exon_dict, overlap_dict, masking)

    composed_paths = _compose_paths(exon_dict, path_dict, number_of_ns)

    # Write
    SeqIO.write(format="fasta", sequences=composed_paths, handle=fasta_out)
