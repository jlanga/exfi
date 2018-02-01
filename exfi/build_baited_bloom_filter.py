#!/usr/bin/env python3

"""
Build a baited bloom filters.

On the one hand, biobloomtools have the functionality of:
    - building bloom filters with biobloommaker
    - classifying reads with respect to such Bloom filters with
biobloomcategorizer

On the other hand, abyss contains tools to produce independent Bloom filters
from biobloomtools with three important differences:
    - levels: can discard low frequency kmers
    - abyss-bloom kmers can report positive kmers in bed format
    - abyss-sealer can fix small (SNPs, INDELS) and not so small gaps (cross
between exons)

We define multiple helper functions:
    - To compose commands:
        - _get_biobloommaker_command
        - _get_categorize_command
        - _get_build_bf_command

And the end product is build_baited_bloom_filter, that stitches together the
pipeline:
1. Build the transcriptomic Bloom filter
2. Get the reads that seem related to the transcriptome
3. Build a bigger bloom filter for downstream analysis with abyss-bloom and
abyss-sealer.
"""

import logging
import shutil
import os

def _get_biobloommaker_command(args, output_dir):
    """Helper function to compose the command to execute"""
    biobloommaker_path = shutil.which('biobloommaker')
    build_transcriptome_bf = [
        biobloommaker_path,
        '--file_prefix', "transcriptome",
        '--output_dir', output_dir,
        '--threads', str(args["threads"]),
        '--kmer_size', str(args["kmer"]),
        args["input_fasta"]
    ]
    return build_transcriptome_bf



def _get_categorize_command(args, output_dir):
    """Helper function to compose the categorize command"""
    categorize_path = shutil.which('biobloomcategorizer')
    categorize = [
        categorize_path,
        '--prefix', output_dir + '/categories',
        '--filter_files', output_dir+ '/transcriptome.bf',
        '--threads', str(args["threads"]),
        '--score', str(args["kmer"]),
        '--fa',
        '--stdout_filter', 'transcriptome',
    ] + args["reads"]
    return categorize



def _get_build_bf_command(args, out_fn):
    """Helper function to compose command to get the final Bloom Filter"""
    abyss_bloom_path = shutil.which('abyss-bloom')
    build_bf = [
        abyss_bloom_path, 'build',
        '--verbose',
        '--kmer', str(args["kmer"]),
        '--bloom-size', args["bloom_size"],
        '--levels', str(args["levels"]),
        '--threads', str(args["threads"]),
        args["output_bloom"]
    ] + out_fn
    return build_bf



def build_baited_bloom_filter(args):
    '''(dict) -> None

    Run the build_baited_bloom_filter pipeline:
    - Build a secondary Bloom filter of the transcriptome with biobloommaker
    - Categorize reads with biobloomcategorizer and pipe to abyss-bloom build
        to build the primary Bloom filter

    Inputs are:
    - A transcriptome in Fasta format
    - A kmer size (better > 25)
    - The total memory to use in the primary Bloom filter
    - The number of levels to use in the Bloom filter (better to use 2)
    - The output file name for the primary Bloom filter
    - The number of threads to use in construction
    - A list of read filenames
    '''

    # Imports
    from subprocess import Popen, PIPE
    from os.path import dirname, abspath

    output_dir = dirname(abspath(args["output_bloom"]))
    # Convert single read library to list
    if isinstance(args["reads"], str):
        args["reads"] = [args["reads"]]

    # Prepare the commands
    build_transcriptome_bf = _get_biobloommaker_command(args, output_dir)
    categorize = _get_categorize_command(args, output_dir)

    build_bf = _get_build_bf_command(args, ["/dev/stdin"])

    # Run the pipeline
    logging.info(
        "\n\nRunning command: %s\n", " ".join(build_transcriptome_bf)
    )

    # Create links to /dev/null for categories_{match,nomatch,multi}.fa and summary
    os.system("ln -s /dev/null " + output_dir + "/categories_transcriptome.fa")
    os.system("ln -s /dev/null " + output_dir + "/categories_noMatch.fa")
    os.system("ln -s /dev/null " + output_dir + "/categories_multiMatch.fa")
    os.system("ln -s /dev/null " + output_dir + "/categories_summary.tsv")

    p_build_transcriptome_bf = Popen(build_transcriptome_bf, shell=False)

    p_build_transcriptome_bf.wait()

    logging.info(
        "\n\nRunning command: %s | %s\n", " ".join(categorize), " ".join(build_bf)
    )

    p_categorize = Popen(categorize, stdout=PIPE, shell=False)
    p_build_bf = Popen(build_bf, stdin=p_categorize.stdout, shell=False)

    p_categorize.stdout.close()
    p_categorize.wait()
    p_build_bf.wait()

    # Clean up files from biobloommaker
    if os.path.isfile(output_dir + "/transcriptome.bf"):
        os.remove(output_dir + "/transcriptome.bf")
    if os.path.isfile(output_dir + "/transcriptome.txt"):
        os.remove(output_dir + "/transcriptome.txt")

    # Clean up files from biobloomcategorizer
    os.remove(output_dir + "/categories_transcriptome.fa")
    os.remove(output_dir + "/categories_noMatch.fa")
    os.remove(output_dir + "/categories_multiMatch.fa")
    os.remove(output_dir + "/categories_summary.tsv")
