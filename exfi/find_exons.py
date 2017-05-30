#!/usr/bin/env python3

# Import everything
from subprocess import Popen, PIPE
from sys import stdin, stdout, stderr
from Bio import SeqIO
from exfi.tab_to_seqrecord import tab_to_seqrecord
from exfi.reduce_exons import reduce_exons


def _abyss_bloom_kmers_command(kmer, bloom_filter_fn, transcriptome_fn):
    """Test all kmers in transcriptome_fn"""
    return [ "abyss-bloom", "kmers",
            "--kmer", str(kmer),
            "--verbose",
            "--bed",
            bloom_filter_fn,
            transcriptome_fn
        ]

def _bedtools_merge_command(kmer):
    """Merge overlapping kmers"""
    return ["bedtools", "merge",
        "-d", str(- kmer + 2)
    ]

def _bedtools_getfasta_command(transcriptome_fn):
    """Get transcriptid:coordinates TAB sequence"""
    return [ "bedtools", "getfasta",
        "-fi", transcriptome_fn,
        "-bed", "-",
        "-tab"
    ]

def _find_exons_pipeline(kmer, bloom_filter_fn, transcriptome_fn):
    # Prepare the commands
    abyss_bloom_kmers = _abyss_bloom_kmers_command(
        kmer, bloom_filter_fn, transcriptome_fn
    )
    bedtools_merge = _bedtools_merge_command(kmer)
    bedtools_getfasta = _bedtools_getfasta_command(transcriptome_fn)


    p1 = Popen(abyss_bloom_kmers, stdout= PIPE)
    p2 = Popen(bedtools_merge, stdin= p1.stdout, stdout= PIPE)
    p3 = Popen(bedtools_getfasta, stdin= p2.stdout, stdout= PIPE)

    # Manage all processes properly
    p1.stdout.close()
    p2.stdout.close()
    pipeline_output_raw = p3.communicate()
    pipeline_output_stdout = pipeline_output_raw[0].decode()
    del pipeline_output_raw
    pipeline_output = pipeline_output_stdout.split("\n")
    p1.wait()
    p2.wait()

    for line in pipeline_output:
        yield line



def find_exons(transcriptome_fn, kmer, bloom_filter_fn, output_fasta):
    """(str, int, str) -> str
    Run the find exons pipeline:
        - abyss-bloom kmers: Test all kmers in the transcriptome
        - bedtools merge: Check overlap and merge
        - bedtools getfasta: convert bed to fasta
    Inputs are:
        - transcriptome_fn: fasta with the transcriptome
        - kmer: int with the kmer length
        - bloom_filter_fn: a bloom filter from abyss-bloom build. Use same
            kmer as above
        - output_fasta: fasta with the different exons
    """

    pipeline_output = _find_exons_pipeline(
        kmer, bloom_filter_fn, transcriptome_fn
    )

    # Process the results from the pipes
    seqrecords = tab_to_seqrecord(pipeline_output)
    exons = reduce_exons(seqrecords)  # Collapse identical exons into one
    SeqIO.write(
        sequences= exons,
        handle= output_fasta,
        format= "fasta"
    )
