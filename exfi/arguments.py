#!/usr/bin/env python3

'''exfi.arguments.py: submodule to deal with the argument parsers'''

import argparse

from exfi import __version__

DEFAULT_K = 25
DEFAULT_BLOOM_SIZE = '500M'
DEFAULT_LEVELS = 2
DEFAULT_THREADS = 1
DEFAULT_MAX_FP_BASES = 5
DEFAULT_MAX_OVERLAP = 10
DEFAULT_MAX_GAP_SIZE = 10
DEFAULT_NUMBER_OF_NS = 100
DEFAULT_COLLAPSE = False

EPILOG = 'Jorge Langa. Send issues and pull requests to github.com/jlanga/exfi'


def add_common_args(parser):
    """Add the common args to a parser: version, verbose, debug"""

    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {version}'.format(
            version=__version__
        )
    )

    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        dest="verbose",
        help="Increase output verbosity"
    )

    parser.add_argument(
        "-d", "--debug",
        action="store_true",
        dest="debug",
        help="Log everything!"
    )


def add_input_fasta(parser):
    """Add the input fasta parser"""
    parser.add_argument(
        '--input-fasta', '-f',
        type=str,
        required=True,
        help='Input transcriptome in FASTA format',
        dest='fasta',
        metavar='FILE'
    )


def add_kmer_size(parser):
    """add the kmer size parser"""
    parser.add_argument(
        '--kmer-size', '-k',
        type=int,
        required=False,
        help=f'The size of the k-mer [{DEFAULT_K}]',
        dest='kmer',
        metavar='INT',
        default=DEFAULT_K
    )


def add_bloom_size(parser):
    """Add the bloom size parser"""
    parser.add_argument(
        '--bloom-size', '-b',
        type=str,
        required=False,
        help=f"Size of the Bloom filter [{DEFAULT_BLOOM_SIZE}]. This is the "
             " total size. The final Bloom filter will be size / levels.",
        dest="bloom_size",
        metavar='STR',
        default=DEFAULT_BLOOM_SIZE
    )


def add_levels(parser):
    """Add the bloom filter levels parser"""
    parser.add_argument(
        '--levels', '-l',
        type=int,
        required=False,
        help=f'Build a cascading bloom filter with N levels and output the '
             f'last level [{DEFAULT_LEVELS}]',
        dest='levels',
        metavar='INT',
        default=DEFAULT_LEVELS
    )


def add_threads(parser):
    """Add the threads parser"""
    parser.add_argument(
        "-t", "--threads",
        type=int,
        help=f'number of threads to use when possible [{DEFAULT_THREADS}]',
        dest="threads",
        metavar="INT",
        default=DEFAULT_THREADS
    )


def add_output_bloom(parser):
    """Add the output bloom parser"""
    parser.add_argument(
        '--output-bloom', '-o',
        type=str,
        required=True,
        help='Path to write the resulting Bloom filter',
        dest="bloom",
        metavar="FILE"
    )


def add_input_bloom(parser):
    """Add the input bloom parser"""
    parser.add_argument(
        '--input-bloom',
        '-b',
        type=str,
        required=True,
        help='Bloom filter with genomic sequences (from '
             'build_baited_bloom_filter or abyss-bloom)',
        dest='bloom',
        metavar='BLOOM'
    )


def add_max_fp_bases(parser):
    """Add the add_max_fp_bases parser"""
    parser.add_argument(
        '--max-fp-bases',
        '-m',
        type=int,
        required=False,
        help='Maximum number of consecutive false positives '
             f'[{DEFAULT_MAX_FP_BASES}]',
        dest="max_fp_bases",
        metavar="INT",
        default=DEFAULT_MAX_FP_BASES
    )


def add_max_overlap(parser):
    """Add the max overlap parser"""
    parser.add_argument(
        '--max-overlap',
        '-l',
        type=int,
        required=False,
        help='Maximum overlap in bp between consecutive exons (0 <= l <= k) '
             f' [{DEFAULT_MAX_OVERLAP}]',
        dest='max_overlap',
        metavar='INT',
        default=DEFAULT_MAX_OVERLAP
    )


def add_max_gap_size(parser):
    """Add the max gap size parser"""
    parser.add_argument(
        '--max-gap-size',
        '-g',
        type=int,
        required=False,
        help="Maximum gap size between predicted exons to try to fill with "
             f"sealer [{DEFAULT_MAX_GAP_SIZE}]",
        dest="max_gap_size",
        metavar='INT',
        default=DEFAULT_MAX_GAP_SIZE
    )


def add_output_gfa(parser):
    '''Add the output gfa parser'''
    parser.add_argument(
        '--output-gfa',
        '-o',
        type=str,
        required=True,
        help='Path to output GFA1 file (the splice graph)',
        dest="gfa1",
        metavar="FILE"
    )


def add_collapse(parser):
    """Add the collapse parser"""
    parser.add_argument(
        '--collapse',
        '-c',
        help='Collapse splice graph by exon sequence [False]',
        dest='collapse',
        action="store_true",
        default=DEFAULT_COLLAPSE
    )


def add_correct(parser):
    """Add the correct parser"""
    parser.add_argument(
        '--correct', '-C',
        help='Correct splice graph by using sealer between exons that seem '
             'nearby [False]',
        dest='correct',
        action="store_true"
    )


def add_polish(parser):
    """Add the polish parser"""
    parser.add_argument(
        "-p", "--polish",
        action="store_true",
        dest="polish",
        help="Polish overlaps in which a AG-GT signal is detected"
    )


def add_input_gfa(parser):
    """Add the input_gfa parser"""
    parser.add_argument(
        '--input-gfa',
        '-i',
        type=str,
        required=True,
        help='Input splice graph in GFA1 format (the results from '
             'build_splicegraph)',
        dest='gfa1',
        metavar='FILE'
    )


def add_output_fasta(parser):
    """Add the output_fasta parser"""
    parser.add_argument(
        '--output-fasta', '-o',
        type=str,
        required=True,
        help='Path to output fasta file',
        dest="fasta",
        metavar="FILE"
    )


def add_soft_mask_overlaps(parser):
    """Add the soft_mask_overlaps parser"""
    parser.add_argument(
        '--soft-mask-overlaps', '-s',
        action='store_true',
        required=False,
        help='Mask overlaps as lowercase nucleotides [False]',
        dest="soft_mask_overlaps",
        default=False
    )


def add_hard_mask_overlaps(parser):
    """"Add the hard_mask_overlaps parser"""
    parser.add_argument(
        '--hard-mask-overlaps', '-m',
        action='store_true',
        required=False,
        help='Mask overlaps as Ns [False]',
        dest="hard_mask_overlaps",
        default=False
    )

def add_to_gapped_transcript(parser):
    """Add the gapped_transcript parser"""
    parser.add_argument(
        '--gapped-transcript', '-g',
        action='store_true',
        required=False,
        help='Convert to gapped transcript instead of separated exons',
        dest='gapped_transcript',
        default=False
    )

def add_number_of_ns(parser):
    """Add the number_of_ns parser"""
    parser.add_argument(
        '--number-of-ns', '-n',
        type=int,
        required=False,
        help=f'Put this number of N between exons [{DEFAULT_NUMBER_OF_NS}]',
        dest="number_of_ns",
        metavar="INT",
        default=DEFAULT_NUMBER_OF_NS
    )



def add_input_splice_graph(parser):
    """Add the input_splice_graph parser"""
    parser.add_argument(
        '--input-splice-graph', '-i',
        type=str,
        required=True,
        help='Input splice graph in BED, GFA or GFF3 format',
        dest='input_splice_graph',
        metavar='FILE'
    )

def add_input_gff3(parser):
    """Add the input_gff3 parser"""
    parser.add_argument(
        '--input-gff3', '-t',
        type=str,
        required=True,
        help='Input splice graph in GFF3 format (a genome assembly annotation)',
        dest='input_gff3',
        metavar='FILE'
    )

def add_gff3_type(parser):
    """Add the gff3_type parser"""
    parser.add_argument(
        '--type-gff3', '-m',
        type=str,
        help='Source of the GFF3 type (ensembl or gmap)',
        choices=['ensembl', 'gmap'],
        default='ensembl',
        metavar='STR',
        dest='gff3_type'
    )

def add_input_fraction(parser):
    """Add the fraction parser"""
    parser.add_argument(
        '--simmilarity-fraction', '-s',
        type=float,
        help='Minimum simmilarity between exons to consider them equal',
        default=0.95,
        metavar='FLOAT',
        dest='fraction'
    )


def build_baited_bloom_filter_args():
    """Create the parser for build_baited_bloom_filter"""
    parser = argparse.ArgumentParser(
        usage='build_baited_bloom_filter '
              '-i transcriptome.fa '
              '-o bloom_filter.bf '
              '-k 30 '
              'reads1.fq ... readsn.fq',
        description=\
            'Build a Bloom filter with reads that have at least one kmer in '
            'the transcriptome.',
        epilog=EPILOG
    )
    add_common_args(parser)
    add_input_fasta(parser)
    add_kmer_size(parser)
    add_bloom_size(parser)
    add_levels(parser)
    add_threads(parser)
    add_output_bloom(parser)
    parser.add_argument(
        metavar='reads',
        type=str,
        nargs='+',
        help='FASTA/Q files (gz or not)',
        dest='reads'
    )

    return parser


def build_splice_graph_args():
    """Create the parser for build_splice_graph"""
    parser = argparse.ArgumentParser(
        usage='build_splice_graph -i transcriptome.fa -b bloom_filter.bf -k 30 '
              '-o exome.gfa',
        description='Store the predicted exome in GFA format',
        epilog=EPILOG
    )
    add_common_args(parser)
    add_input_fasta(parser)
    add_input_bloom(parser)
    add_kmer_size(parser)
    add_max_fp_bases(parser)
    add_max_overlap(parser)
    add_max_gap_size(parser)
    add_output_gfa(parser)
    add_correct(parser)
    add_polish(parser)
    add_threads(parser)
    add_collapse(parser)
    return parser


def gfa1_to_fasta_args():
    """Create the gfa1_to_exons parser"""
    parser = argparse.ArgumentParser(
        usage='gfa_to_fasta -i splice_graph.gfa -o exons.fa',
        description='Extract the exons from a splice graph in GFA format. '
                    'Optionally mask overlaps between consecutive exons',
        epilog=EPILOG
    )
    add_common_args(parser)
    add_input_gfa(parser)
    add_output_fasta(parser)
    add_to_gapped_transcript(parser)
    add_number_of_ns(parser)
    add_soft_mask_overlaps(parser)
    add_hard_mask_overlaps(parser)
    return parser



def compare_to_gff_args():
    """Create the compare_gfa_to_gff parser"""
    parser = argparse.ArgumentParser(
        usage='compare_gff_to_gfa -i splice_graph.gfa -t annotation.gff3 -m '
              'ensembl -f 0.95',
        description=\
            'Compare a splice graph in GFA1/BED/GFF3 format against a genome annotation '
            'in gff3 format.',
        epilog=EPILOG
    )
    add_common_args(parser)
    add_input_splice_graph(parser)
    add_input_gff3(parser)
    add_gff3_type(parser)
    add_input_fraction(parser)

    return parser
