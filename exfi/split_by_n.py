#!/usr/bin/env python3

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def split_by_n(iterable_seqs):
    """(generator)-> generator 
    Split a record by the Ns it contains into multiple Seqs
    Each new seq.id  is seq.id::subNUMBER
    """
    # dict to know
    contig2sequence = {}
    for record in iterable_seqs:
        subrecord_id = 0
        record_id = record.id

        splitted_sequences = [x for x in str(record.seq).split("N") if x != '']
        
        for sequence in splitted_sequences:
            subrecord_id += 1
            subrecord = SeqRecord(
                Seq(
                    sequence,
                    IUPAC.IUPACAmbiguousDNA
                ),
                id = record_id + "_s" + str(subrecord_id), 
                name = "",
                description = " "
            )

            yield subrecord