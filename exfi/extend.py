#!/usr/bin/env python3

def extend_left(records, kmer):
    """(iterable of SeqRecord, int) -> generator

    Given a Seq and an integer representing the length of the kmer,
    yield the four possible extensions to the left of the record.
    Records shorter than kmer - 1 won't be extended or returned.
    """
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    assert isinstance(kmer, int) and kmer > 0, \
        "Incompatible kmer {}".format(kmer)

    records = (record for record in records if len(record) >= kmer -1)

    for record in records:
        for letter in "ACGT":
            new_record = SeqRecord(
                id = record.id,
                seq = record.seq,
                description = ""
            )
            new_record.id += "l" + letter
            new_record.seq = Seq(letter + str(new_record.seq)[:kmer-1])
            yield new_record


def extend_right(records, kmer):
    """(iterable of Seqrecords, int) -> generator of SeqRecords

    Given a Seq and an integer representing the length of the kmer,
    yield the four possible extensions to the right of the record.
    Records shorter than kmer - 1 won't be extended or returned.
    """
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    assert isinstance(kmer, int) and kmer > 0, \
        "Incompatible kmer {}".format(kmer)

    records = (record for record in records if len(record) >= kmer -1)

    for record in records:
        for letter in "ACGT":
            new_record = SeqRecord(
                id = record.id,
                seq = record.seq,
                description = ""
            )
            new_record.id += "r" + letter
            new_record.seq = Seq(str(new_record.seq)[-kmer+1:] + letter)
            yield new_record
