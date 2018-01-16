#!/usr/bin/env python3

"""
Auxiliary functions and classes for testing
"""

import tempfile
import shutil
from subprocess import Popen, PIPE

import networkx as nx
from Bio import SeqIO

from exfi.find_exons import \
    _process_output, \
    _get_fasta, \
    _find_exons_pipeline

from exfi.build_baited_bloom_filter import \
    _get_build_bf_command






def _command_to_list(command):
    """Execute command and return output as list of strings"""
    process = Popen(command, stdout=PIPE, shell=False)
    results = list(_process_output(process))
    return results

def _fasta_to_dict(filename):
    """SeqIO.index wrapper for fasta files"""
    return SeqIO.index(filename=filename, format="fasta")


def _fasta_to_list(filename):
    """SeqIO.parse wrapper for fasta files"""
    return list(SeqIO.parse(handle=filename, format="fasta"))

def _getfasta_to_list(transcriptome_dict, iterable_of_bed):
    """Convert to a list the generator from getfasta"""
    return list(_get_fasta(transcriptome_dict, iterable_of_bed))

def _silent_popen(command):
    """Create a Popen with no stderr and stdout"""
    return Popen(
        command,
        stdout=open("/dev/null", 'w'),
        stderr=open("/dev/null", 'w'),
        shell=False
    )

def _bf_and_process(reads_fns, transcriptome_fn):
    """(list of str, str) -> list

    Build the BF and process the reads
    """
    tmp_dir = tempfile.mkdtemp()
    tmp_bf = tmp_dir + "/transcriptome_noreads.bf"
    args = {
        "kmer": 30,
        "bloom_size": "100M",
        "levels": 1,
        "threads": 1,
        "input_bloom": tmp_bf,
        "output_bloom": tmp_bf,
        "reads": reads_fns,
        "input_fasta": transcriptome_fn,
        "max_fp_bases": 5,
        "max_overlap": 10
    }
    command = _get_build_bf_command(args, reads_fns)
    process = _silent_popen(command)
    process.wait()
    results = _find_exons_pipeline(args)
    shutil.rmtree(tmp_dir)
    return list(results)


class CustomAssertions:
    """
    Custom assertions not covered in unittest:
    - assertEqualListOfSeqrecords
    """
    @classmethod
    def assertEqualListOfSeqrecords(self, records1, records2):
        """
        Check if each element of list_of_seqrecords1 is exactly equal to each one of
        list_of_seqrecords2.
        """
        # pylint: disable=invalid-name, bad-classmethod-argument
        length_1 = len(records1)
        length_2 = len(records2)
        if length_1 != length_2:
            raise AssertionError(
                'Lengths differ:\n {len_1} != {len_2}'.format(
                    len_1=length_1,
                    len_2=length_2
                )
            )
        else:
            for i in range(length_1):
                record1 = records1[i]
                record2 = records2[i]
                if record1.id != record2.id or record1.seq != record2.seq:
                    raise AssertionError(
                        'Records at position {i} differ:\n{id1} : {seq1}\n{id2} : {seq2}'.format(
                            i=i,
                            id1=record1.id,
                            seq1=record1.seq,
                            id2=record2.id,
                            seq2=record2.seq
                        )
                    )
            return True

    @classmethod
    def assertEqualSpliceGraphs(self, sg1, sg2):
        """
        Check if two splice graph are equal:
        - are isomorphic
        - same coordinates
        - same overlaps
        """
        # pylint: disable=invalid-name,bad-classmethod-argument
        coordinates1 = nx.get_node_attributes(G=sg1, name="coordinates")
        coordinates2 = nx.get_node_attributes(G=sg1, name="coordinates")
        overlaps1 = nx.get_edge_attributes(G=sg1, name="overlaps")
        overlaps2 = nx.get_edge_attributes(G=sg1, name="overlaps")
        same_coordinates = coordinates1 == coordinates2
        same_overlaps = overlaps1 == overlaps2
        are_isomorphic = nx.is_isomorphic(sg1, sg2)
        if are_isomorphic and same_coordinates and same_overlaps:
            return True
        return False
