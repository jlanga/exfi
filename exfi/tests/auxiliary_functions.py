#!/usr/bin/env python3

"""
Auxiliary functions and classes for testing
"""

from exfi.find_exons import \
    _process_output, \
    _get_fasta, \
    _find_exons_pipeline

from exfi.build_baited_bloom_filter import \
    _get_build_bf_command

from subprocess import Popen, PIPE
from Bio import SeqIO

import tempfile
import shutil

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
    return Popen(command,
        stdout=open("/dev/null", 'w'),
        stderr=open("/dev/null", 'w'),
        shell=False
    )

def _bf_and_process(reads_fns, transcriptome_fn):
    """Build the BF and process the reads"""
    tmp_dir = tempfile.mkdtemp()
    tmp_bf = tmp_dir + "/transcriptome_noreads.bf"
    command = _get_build_bf_command("30", "100M", "1", "1", tmp_bf, reads_fns)
    process = _silent_popen(command)
    process.wait()
    results = _find_exons_pipeline(
        kmer=30,
        bloom_filter_fn=tmp_bf,
        transcriptome_fn=transcriptome_fn,
        max_fp_bases=5
    )
    shutil.rmtree(tmp_dir)
    return list(results)


class CustomAssertions:

    @classmethod
    def assertEqualListOfSeqrecords(self, records1, records2):
        """
        Check if each element of list_of_seqrecords1 is exactly equal to each one of
        list_of_seqrecords2.
        """
        n1 = len(records1)
        n2 = len(records2)
        if n1 != n2:
            raise AssertionError(
                'Lengths differ:\n {len_1} != {len_2}'.format(
                    len_1 = n1,
                    len_2 = n2
                )
            )
        else:
            for i in range(n1):
                record1=records1[i]
                record2=records2[i]
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
