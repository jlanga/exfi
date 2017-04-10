#!/usr/bin/env python3

def reduce_exons(iterable_of_seqrecords):
    """(iterable of SeqRecords) -> iterable of SeqRecords
    
    Reduce redundacy of the iterable seqrecords by removing repeated exons
    and storing graph info into the description header
    """

    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    # Make an association between the exon sequence and info from the transcriptome
    # keys are sequences (as str) and values are 
    # [new_exon_id, old_exon1_id, ..., old_exonN_id]
    seq_to_exon = dict()

    for record in iterable_of_seqrecords:

        # Split record into relevant info
        sequence = str(record.seq)
    
        # Initialize the sequence
        if str(record.seq) not in seq_to_exon:
            new_exon_id = 'EXON{:011d}'.format(len(seq_to_exon)+1)
            seq_to_exon[sequence] = [new_exon_id]

        # Just add the other info
        seq_to_exon[sequence].append(record.id)

    # revert the dictionary and return sequences
    # Dict has the structure exon_id : seqrecord
    
    exonid_to_seqrecord = {}

    for sequence in seq_to_exon:

        info = seq_to_exon[sequence]

        identifier = info[0]
        description = " ".join(info[1:])

        record = SeqRecord(
            id = identifier,
            description = description,
            seq = Seq(sequence)
        )

        exonid_to_seqrecord[identifier] = record

    # Return results ordered by exon_id
    for identifier in sorted(exonid_to_seqrecord):
        yield exonid_to_seqrecord[identifier]