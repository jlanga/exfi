#!/usr/bin/env python3

"""exfi.io.gfa1_to_exons.py: submodule to read a gfa1, extract the exons and
store it in fasta format"""

import logging

from exfi.io.read_gfa import read_gfa1
from exfi.io.masking import mask, cigar_to_int

def gfa1_to_fasta(
        fasta_out, gfa1_in, masking='none', transcripts=False,
        number_of_ns=100):
    """Extract the exons from a GFA1 file and store them in fasta format

    if gapped_transcripts is True, the result is a gapped transcript, where the
    exons are separated by number_of_ns Ns.
    """
    logging.info('Converting GFA1 to fasta')

    with open(fasta_out, 'w') as fasta:

        logging.info('Reading the GFA1 file')
        gfa1 = read_gfa1(gfa1_in)

        logging.info('Computing node2sequence')
        node2sequence = gfa1['segments']\
            .drop(columns='record_type')

        if node2sequence.shape[0] == 0:
            return

        logging.info('Computing edge2overlap')
        edge2overlap = gfa1['links']\
            .drop(columns=['record_type', 'from_orient', 'to_orient'])\
            .rename(columns={
                'from': 'u', 'to': 'v', 'overlap': 'overlap_cigar'
            })
        logging.info('Computing overlap from CIGAR to int')
        edge2overlap['overlap'] = edge2overlap.overlap_cigar.map(cigar_to_int)

        logging.info('Masking (if necessary)')
        node2sequence = mask(
            node2sequence=node2sequence,
            edge2overlap=edge2overlap,
            masking=masking
        )

        if transcripts:
            separator = number_of_ns * 'N'

            logging.info('Computing path2nodes')
            path2nodes = gfa1['paths']\
                .drop(columns=['record_type', 'overlaps'])
            path2nodes.segment_names = path2nodes\
                .segment_names.str.replace('+', '')

            logging.info('Converting node2sequence to dict')
            node2sequence_dict = node2sequence\
                .set_index('name')\
                .to_dict()['sequence']

            logging.info('Composing the gapped sequence')
            path2nodes["gapped_sequence"] = path2nodes\
                .segment_names\
                .str.split(',')\
                .map(
                    lambda x: separator.join([node2sequence_dict[y] for y in x])
                )

            logging.info('Composing fasta sequences')
            path2nodes["fasta"] = \
                ">" + path2nodes.path_name + " " + path2nodes.segment_names + \
                "\n" + \
                path2nodes.gapped_sequence

            logging.info('Dumping fasta to disk')
            path2nodes.fasta.values.tofile(fasta, sep="\n", format="%s")
            fasta.write("\n")  # Final end line
            logging.info('Done')

        else:
            logging.info('Composing fasta sequences')
            node2sequence["fasta"] = \
                ">" + node2sequence["name"] + "\n" + \
                node2sequence["sequence"]

            logging.info('Dumping fasta to disk')
            node2sequence.fasta.values.tofile(fasta, sep="\n", format="%s")
            fasta.write("\n")  # Final end line
            logging.info('Done')
