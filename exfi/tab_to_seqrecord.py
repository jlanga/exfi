#!/usr/bin/env python3


from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def tab_to_seqrecord(iterable):
    for element in iterable:
        if element:
            identifier, sequence = element.strip().split("\t")
            new_identifier = identifier.rsplit(":")[0]
            description = identifier
            yield SeqRecord(
                id=new_identifier,
                seq=Seq(sequence),
                description=description
            )
