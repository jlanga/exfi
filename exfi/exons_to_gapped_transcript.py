#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def build_transcript_to_exon_dict(exons):
    """
    (Iterable of SeqRecords) -> dict

    Build a dictionary where the keys are the transcript_id, and the values
    a tuple of exon_ids, sorted by their coordinates
    """
    # Build a dict of transcript_id : (exon_id1, ..., exon_idN)
    # Exons have as identifier exon_id and trx_id:start-end as description
    transcript_to_exons = {}

    for exon in exons:

        exon_id = exon.id

        # Get the string of descriptions
        coordinates = exon.description.split(" ")[1:]

        # Process every substring, getting the transcript_id and the start of
        # the exon wrt the transcript
        for coordinate in coordinates:

            # Get tr_id and start from tr_id:start-end
            transcript_id = coordinate.split(":")[0]
            start = int(coordinate.rsplit(":")[1].split("-")[0])

            # Initialize value if not there
            if transcript_id not in transcript_to_exons:
                transcript_to_exons[transcript_id] = []

            # Append
            transcript_to_exons[transcript_id].append(
                (exon_id, start)
            )

    #Sort the transcript_to_exons values by start position, and delete it
    for transcript_id in transcript_to_exons:
        # Get the elements in such transcript
        exons_in_transcript = transcript_to_exons[transcript_id]
        exons_in_transcript = sorted( # Sort them
            exons_in_transcript,
            key = lambda x: x[1]
        )
        exons_in_transcript = [x[0] for x in exons_in_transcript] # Get only the exon_id and dump the rank
        transcript_to_exons[transcript_id] = tuple(exons_in_transcript)

    return transcript_to_exons


def exon_dict_to_gapped_transcript(transcript_to_exons, exome_fn, number_of_ns=100):
    """(dict of SeqRecords) -> iterable of SeqRecords

    Paste together all exons of a transcript, separated by Ns.
    """
    # Index so we can get subsequences (exons) from it
    exome = SeqIO.index(filename=exome_fn, format="fasta")
    ns = "N" * number_of_ns

    for transcript_id in transcript_to_exons.keys():
        # Get all the exon_ids required that compose the transcript
        exon_ids = list(element for element in transcript_to_exons[transcript_id])
        # Get all sequences as
        sequences = [str(exome[exon_id].seq) for exon_id in exon_ids if exon_id in exome]
        sequence = ns.join(sequences)
        record = SeqRecord(
            id=transcript_id,
            seq=Seq(sequence),
            description="",
            name=""
        )
        yield record
