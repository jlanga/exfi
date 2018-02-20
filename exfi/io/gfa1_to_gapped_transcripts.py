#!/usr/bin/env python3

"""exfi.io.gfa1_to_gapped_transcript.py: submodule to convert a GFA1 file to a fasta file where the
spaces between predicted exons are filled with a string of Ns.
"""

import logging

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from exfi.io.masking import _mask
from exfi.io.read_gfa1 import read_gfa1



def _compose_paths(exon_dict: dict, path_dict: dict, number_of_ns: int = 100) -> SeqRecord:
    """Compose and return each gapped transcript.

    :param exon_dict: dict of exons: {exon_id: ((seq1, start1, end1,), ...)}.
    :param path_dict: dict of paths: {transcript1: (exon1, ..., exonN)}.
    :param number_of_ns: number of Ns to write between each exon.
    """
    logging.info("\tComposing paths")
    chunk_of_ns = "N" * number_of_ns
    for transcript_id, exon_list in sorted(path_dict.items()):
        exon_seqs = [str(exon_dict[exon_id]) for exon_id in exon_list]
        yield SeqRecord(
            id=transcript_id,
            description=",".join(exon_list),
            seq=Seq(chunk_of_ns.join(exon_seqs))
        )


def gfa1_to_gapped_transcripts(
        gfa_in: str, fasta_out: str, number_of_ns: int = 100, masking: str = "none") -> None:
    """Write gapped transcripts as fasta from GFA1 file

    :param str gfa_in: Path to input GFA1 file.
    :param str fasta_out: Patho to output FASTA file.
    :param int number_of_ns: Number of Ns to be written between each exon (Default value = 100).
    :param str masking: Type of masking to be applied: none, soft, hard (Default value = "none").
    """

    logging.info("Converting GFA1 file %s to gapped transcript fasta %s", gfa_in, fasta_out)

    # Process
    gfa1 = read_gfa1(gfa_in)
    exon_dict = gfa1["segments"]
    overlap_dict = gfa1["links"]
    path_dict = gfa1["paths"]


    # Mask if necessary
    exon_dict = _mask(exon_dict, overlap_dict, masking)

    composed_paths = _compose_paths(exon_dict, path_dict, number_of_ns)

    # Write
    logging.info("\tWriting to file")
    SeqIO.write(format="fasta", sequences=composed_paths, handle=fasta_out)
