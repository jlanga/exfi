#!/usr/bin/env python3

def exons_to_gapped_transcript(exons, gapped_transcripts):
    '''(iterable_of_seqrecords) -> iterable_of_seqrecords

    Take an iterable of exons and rebuild the constituent gapped transcript
    '''

    for exon in exons:

        exon_id = exon.id
        locations = [x.exon.description.split(" ")[1:]
        