#!/usr/bin/env python3


# Imports
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def _make_seq_to_exon(iterable_of_seqrecords):
    # 1. Make an association between the exon sequence and info from
    # the transcriptome
    # keys are sequences (as str) and values are
    # [new_exon_id, old_exon1_id, ..., old_exonN_id]
    seq_to_exon = dict()  # sequence : [exons_ids with such sequence]

    for record in iterable_of_seqrecords:

        # Split record into relevant info
        sequence = str(record.seq)

        # Initialize the sequence
        if str(record.seq) not in seq_to_exon:
            new_exon_id = 'EXON{:011d}'.format(len(seq_to_exon) + 1)
            seq_to_exon[sequence] = [new_exon_id]

        # Just add the other info
        seq_to_exon[sequence].append(record.description)
    return seq_to_exon


def _exonid_to_seqrecord(seq_to_exon):
    # 2. revert the dictionary and return sequences
    # Dict has the structure exon_id : seqrecord

    exonid_to_seqrecord = {}

    for sequence in seq_to_exon:

        info = seq_to_exon[sequence]

        identifier = info[0]
        description = " ".join(info[1:])

        record = SeqRecord(
            id=identifier,
            description=description,
            seq=Seq(sequence)
        )

        exonid_to_seqrecord[identifier] = record
    return exonid_to_seqrecord


def reduce_exons(iterable_of_seqrecords):
    """(iterable of SeqRecords) -> iterable of SeqRecords

    Reduce redundacy of the iterable seqrecords by removing repeated exons
    and storing graph info into the description header
    """

    seq_to_exon = _make_seq_to_exon(iterable_of_seqrecords)
    exonid_to_seqrecord = _exonid_to_seqrecord(seq_to_exon)

    # Return results ordered by exon_id
    for identifier in sorted(exonid_to_seqrecord):
        yield exonid_to_seqrecord[identifier]
