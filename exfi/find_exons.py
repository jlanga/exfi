#!/usr/bin/env python3

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

    # Import everything
    from subprocess import Popen, PIPE
    from sys import stdin, stdout, stderr
    from Bio import SeqIO
    from exfi.tab_to_seqrecord import tab_to_seqrecord
    from exfi.reduce_exons import reduce_exons

    # Prepare the commands
    abyss_bloom_kmers = [  # Run abyss-bloom kmers
    "abyss-bloom", "kmers",
        "--kmer", str(kmer),
        "--verbose",
        "--bed",
        bloom_filter_fn,
        transcriptome_fn
    ]   

    bedtools_merge = [  # Merge overlapping kmers
    "bedtools", "merge",
        "-d", str(- kmer + 2)
    ]

    bedtools_getfasta = [  # Get transcriptid:coordinates TAB sequence
        "bedtools", "getfasta",
            "-fi", transcriptome_fn,
            "-bed", "-",
            "-tab"
    ]
    
    # Run the pipeline
    # Get all kmers from the transcriptome that are in the genome
    process_abyss_bloom_kmers = Popen(
        abyss_bloom_kmers,
        stdout= PIPE,
    )
    
    # Merge them
    process_bedtools_merge = Popen(
        bedtools_merge,
        stdin= process_abyss_bloom_kmers.stdout,
        stdout= PIPE
    )

    # Build a fasta
    process_bedtools_getfasta = Popen(
        bedtools_getfasta,
        stdin= process_bedtools_merge.stdout,
        stdout= PIPE,
    )
    
    # Manage all processes properly    
    process_abyss_bloom_kmers.stdout.close()
    process_bedtools_merge.stdout.close()
    pipeline_output = process_bedtools_getfasta.communicate()[0].decode().split("\n")
    process_abyss_bloom_kmers.wait()
    process_bedtools_merge.wait()

    # Process the results from the pipes
    seqrecords = tab_to_seqrecord(pipeline_output)
    exons = reduce_exons(seqrecords)
    SeqIO.write(
        sequences= exons,
        handle= output_fasta,
        format= "fasta"
    )