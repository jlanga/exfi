#!/usr/bin/env python3

"""Build a baited bloom filters.

On the one hand, biobloomtools have the functionality of:
    - building bloom filters with biobloommaker
    - classifying reads with respect to such Bloom filters with biobloomcategorizer

On the other hand, abyss contains tools to produce independent Bloom filters from biobloomtools
with three important differences:
    - levels: can discard low frequency kmers
    - abyss-bloom kmers can report positive kmers in bed format
    - abyss-sealer can fix small (SNPs, INDELS) and not so small gaps (cross between exons)

We define multiple helper functions:
    - To compose commands:
        - _get_biobloommaker_command
        - _get_categorize_command
        - _get_build_bf_command

And the end product is build_baited_bloom_filter, that stitches together the pipeline:

1. Build the transcriptomic Bloom filter
2. Get the reads that seem related to the transcriptome
3. Build a bigger bloom filter for downstream analysis with abyss-bloom and abyss-sealer.
"""


import logging
import shutil
import os


def _get_biobloommaker_command(args: dict, output_dir: str) -> list:
    """Helper function to compose the command to execute

    :param dict args: Dict of arguments.
    :param str: output_dir: Output folder.

    :ivar dict args: Dict of arguments. The ones used are {"threads": int, "kmer": int}, where
        - threads is the number of threads to be used
        - kmer is the kmer size
    :ivar str output_dir: Output directory. Must exist and won't be created.


    .. seealso: Where it is used:
    :py:meth:`build_baited_bloom_filter`

    """
    biobloommaker_path = shutil.which('biobloommaker')
    build_transcriptome_bf = [
        biobloommaker_path,
        '--file_prefix', "transcriptome",
        '--output_dir', output_dir,
        '--threads', str(args["threads"]),
        '--kmer_size', str(args["kmer"]),
        args["fasta"]
    ]
    return build_transcriptome_bf


def _get_categorize_command(args: dict, output_dir: str) -> list:
    """Helper function to compose the categorize command

    :arg dict args: Dict of arguments
    :arg str output_dir: The output directory

    :ivar dict args: Dict of arguments. The ones used are:
        - "threads": int: number of threads to be used
        - "kmer": int: the kmer size,
        - "reads": list: list of str of paths to the reads to be categorized
    :ivar str output_dir: Must already exist. Won't be created.

    .. seealso: Where it is used:
        :py:meth: `build_baited_bloom_filter`

    :param dict args: Dict of arguments.
    :param str output_dir: Output directory.

    """
    categorize_path = shutil.which('biobloomcategorizer')
    categorize = [
        categorize_path,
        '--prefix', output_dir + '/categories',
        '--filter_files', output_dir + '/transcriptome.bf',
        '--threads', str(args["threads"]),
        '--score', str(args["kmer"]),
        '--fa',
        '--stdout_filter', 'transcriptome',
    ] + args["reads"]
    return categorize


def _get_build_bf_command(args: dict, in_fn: str) -> list:
    """Helper function to compose command to get the final Bloom Filter

    :arg dict args: Dict of arguments.
    :arg str in_fn: Path to file where the reads will be read

    :ivar dict args: Dict of arguments. The ones used are:
        - "kmer" (int): size of the kmers.
        - "threads" (int): number of threads to be used.
        - "bloom_size" (str): size of the Bloom filter in bytes. K/M/G units can be used.
        - "levels" (int): number of Bloom filters used.
        - "output_bloom" (str): path where to store the Bloom Filter
    :ivar str in_fn: Path to file where the reads will be read. In later methods /dev/stdin is used.



    :param dict args: Dict of arguments.
    :param str in_fn: Input filename.

    ..seealso: Where it is used:
        :py:meth: `build_baited_bloom_filter`

    """
    abyss_bloom_path = shutil.which('abyss-bloom')
    build_bf = [
        abyss_bloom_path, 'build',
        '--verbose',
        '--kmer', str(args["kmer"]),
        '--bloom-size', args["bloom_size"],
        '--levels', str(args["levels"]),
        '--threads', str(args["threads"]),
        args["bloom"]
    ] + in_fn
    return build_bf


def _create_links(output_dir: str) -> None:
    """Create soft links from output_dir/categries-* to /dev/null

    :arg str output_dir: path where data will be stored

    :ivar str output_dir: path where data will be stored. Must exist, won't be created.

    ..seealso: Where it is used:
        :py:meth: `build_baited_bloom_filter`

    :param str output_dir: Output directory.
    """
    os.symlink(os.devnull, output_dir + "/categories_transcriptome.fa")
    os.symlink(os.devnull, output_dir + "/categories_noMatch.fa")
    os.symlink(os.devnull, output_dir + "/categories_multiMatch.fa")
    os.symlink(os.devnull, output_dir + "/categories_summary.tsv")


def _destroy_links(output_dir: str) -> None:
    """Destroy the links created by _create_links

    :arg str output_dir: path where data is being stored.

    :ivar str output_dir: path where data will be stored. Must exist, won't be created.

    ..seealso: Where it is used:
        :py:meth: `build_baited_bloom_filter`

    :param str output_dir: Output directory.
    """
    os.remove(output_dir + "/categories_transcriptome.fa")
    os.remove(output_dir + "/categories_noMatch.fa")
    os.remove(output_dir + "/categories_multiMatch.fa")
    os.remove(output_dir + "/categories_summary.tsv")


def build_baited_bloom_filter(args: dict) -> None:
    """Run the build_baited_bloom_filter pipeline.

    The pipeline works as follows:

        - Build a secondary Bloom filter of the transcriptome with `biobloommaker`.
        - Categorize reads with `biobloomcategorizer` and pipe to `abyss-bloom build` to build the
primary Bloom filter

    :param dict args: Dict of arguments.

    .. note:: Parameter `args` must have the following keys:

        - fasta (str): path to transcriptome in FASTA format.
        - kmer (int): size of the kmers.
        - bloom_size (str):  total size of the Bloom filter(s) in bytes. K/M/G units can be used.
        - levels (int): number of Bloom filters used to count.
        - bloom (str): path where to store the final Bloom filter.
        - threads (int): number of threads to be used.
        - reads (list of str): paths to the reads in FASTA/Q format, gzipped or not.

    .. seealso: Functions used:
        :py:meth: `_get_biobloommaker_command`, `_get_categorize_command`, `_get_build_bf_command`,
        `_create_links`, `_destroy_links`.

    """

    # Imports
    from subprocess import Popen, PIPE
    from os.path import dirname, abspath

    output_dir = dirname(abspath(args["bloom"]))
    # Convert single read library to list
    if isinstance(args["reads"], str):
        args["reads"] = [args["reads"]]

    # Prepare the commands
    build_transcriptome_bf = _get_biobloommaker_command(args, output_dir)
    categorize = _get_categorize_command(args, output_dir)
    build_bf = _get_build_bf_command(args, ["/dev/stdin"])

    # Run the pipeline
    logging.info("\n\nRunning command: %s\n", " ".join(build_transcriptome_bf))

    # Create links to /dev/null for categories_{match,nomatch,multi}.fa and summary
    _create_links(output_dir)

    # Put the porcesses to work together
    p_build_transcriptome_bf = Popen(build_transcriptome_bf, shell=False)
    p_build_transcriptome_bf.wait()
    logging.info("\n\nRunning command: %s | %s\n",
                 " ".join(categorize), " ".join(build_bf))
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
    _destroy_links(output_dir)
