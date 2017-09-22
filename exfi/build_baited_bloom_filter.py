#!/usr/bin/env python3

import shutil

def _get_biobloommaker_command(output_dir, threads, kmer, transcriptome):
    """Helper function to compose the command to execute"""
    biobloommaker_path = shutil.which('biobloommaker')
    build_transcriptome_bf = [
        biobloommaker_path,
        '--file_prefix', "transcriptome",
        '--output_dir', output_dir,
        '--threads', str(threads),
        '--kmer_size', str(kmer),
        transcriptome
    ]
    return build_transcriptome_bf

def _get_categorize_command(output_dir, threads, kmer, reads):
    """Helper function to compose the categorize command"""
    categorize_path = shutil.which('biobloomcategorizer')
    categorize = [
        categorize_path,
        '--prefix', output_dir + '/categories',
        '--filter_files', output_dir + '/transcriptome.bf',
        '--threads', str(threads),
        '--score', str(kmer),
        '--fa',
        '--stdout_filter', 'transcriptome',
    ] + reads
    return categorize

def _get_build_bf_command(kmer, bloom_size, levels, threads, output_bloom, input_fns):
    """Helper function to compose command to get the final Bloom Filter"""
    abyss_bloom_path = shutil.which('abyss-bloom')
    build_bf = [
        abyss_bloom_path, 'build',
        '--verbose',
        '--kmer', str(kmer),
        '--bloom-size', bloom_size,
        '--levels', str(levels),
        '--threads', str(threads),
        output_bloom,
    ] + input_fns
    return build_bf

def build_baited_bloom_filter(
    transcriptome,
    kmer,
    bloom_size,
    levels,
    output_bloom,
    threads,
    reads
):
    '''(str, int, str, int, str, int , list of str) -> None

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
    from sys import stderr
    from os.path import dirname

    output_dir = dirname(output_bloom)
    if output_dir == "": output_dir = "./"

    # Convert single read library to list
    if isinstance(reads, str):
        reads = [reads]

    # Prepare the commands
    build_transcriptome_bf = _get_biobloommaker_command(
        output_dir, threads, kmer, transcriptome
    )

    categorize = _get_categorize_command(
        output_dir, threads, kmer, reads
    )

    build_bf =_get_build_bf_command(
        kmer, bloom_size, levels, threads, output_bloom, ["/dev/stdin"]
    )

    # Run the pipeline
    stderr.write(
        "\n\nRunning command: {command}\n".format(
            command=" ".join(build_transcriptome_bf)
        )
    )

    p_build_transcriptome_bf = Popen(build_transcriptome_bf, shell=False)

    p_build_transcriptome_bf.wait()

    stderr.write(
        "\n\nRunning command: {command1} | {command2}\n".format(
            command1=" ".join(categorize),
            command2=" ".join(build_bf)
        )
    )

    p_categorize = Popen(categorize, stdout=PIPE, shell=False)
    p_build_bf = Popen(build_bf, stdin=p_categorize.stdout, shell=False)

    p_categorize.stdout.close()
    p_categorize.wait()
    p_build_bf.wait()
