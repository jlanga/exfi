#!/usr/bin/env python3

"""
Module to compute positive exons in the bloom filter as follos:
- abyss-bloom kmers to get bed coordinates of each transcriptomic kmer in the BF
- bedtools merge to join overlapping consecutive intervals by k-1bases
- awk to throw away too short intervales (to avoid FPs)
- bedtools merge to join overlapping consecutive intervals by k-max_overlap bases
"""


# Import everything

from typing import \
    Iterable, \
    Tuple

import logging

from subprocess import Popen, PIPE


def _process_output(process: Popen) -> Iterable[Tuple[str, int, int]]:
    """Get lines in bed format from the output of a Popen.

    :param Popen process: Popen object.
    """
    for stdout_line in iter(process.stdout.readline, b''):
        chromosome, start, end = stdout_line.decode().strip().split()
        yield (chromosome, int(start), int(end))
    process.stdout.close()
    process.wait()


def _get_fasta(transcriptome_dict: dict, iterable_of_bed: list) -> Iterable[Tuple[str, str]]:
    """Extract subsequences in trancriptome_fn according to locis.

    :param dict transcriptome_dict: FastaDict of the transcriptome
    :param iterable iterable_of_bed: iterable of Bed3Records.
    """
    for bed in iterable_of_bed:
        chromosome, start, end = bed
        if chromosome in transcriptome_dict:
            seq = transcriptome_dict[chromosome][start:end]
            identifier = "{0}:{1}-{2}".format(chromosome, start, end)
            yield (identifier, seq)


def _find_exons_pipeline(args: dict) -> Iterable[Tuple[str, int, int]]:
    """Find exons according to the Bloom filter -> BED

    Main pipeline:
    - Check every kmer,
    - merge if there is an overlap of k-1 bases,
    - Throw away too short exons (k + max_fp_bases; FP of the BF),
    - merge consecutive exons if they have an ovelap of max_fp_bases

    args = {
        "kmer": int,
        "bloom": str,
        "fasta": str,
        "max_overlap": int,
        "max_fp_bases": int
    }

    :param dict args: arguments for abyss-bloom kmers, bedtools merge and awk

    """
    logging.info("Running the find exons pipeline")
    # Prepare the commands
    c_kmers = [
        "abyss-bloom", "kmers", "--kmer", str(
            args["kmer"]), "--verbose", "--bed",
        args["bloom"], args["fasta"]
    ]
    c_merge1 = ["bedtools", "merge", "-d", str(-args["kmer"] + 2)]
    c_filter = ["awk", "$3 - $2 >= {min_length}".format(
        min_length=args["kmer"] + args["max_fp_bases"]
    )]
    c_merge2 = ["bedtools", "merge", "-d", str(-args["max_overlap"])]
    # Run all of them streamlined
    p_kmers = Popen(c_kmers, stdout=PIPE)
    p_merge1 = Popen(c_merge1, stdin=p_kmers.stdout, stdout=PIPE)
    p_filter = Popen(c_filter, stdin=p_merge1.stdout, stdout=PIPE)
    p_merge2 = Popen(c_merge2, stdin=p_filter.stdout, stdout=PIPE)
    p_kmers.stdout.close()
    yield from _process_output(p_merge2)
