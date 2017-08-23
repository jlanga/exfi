#!/usr/bin/env python3


# Import everything
from subprocess import Popen, PIPE
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from exfi.reduce_exons import reduce_exons


def _process_output(process):
    """Get lines from the output of a Popen."""
    for stdout_line in iter(process.stdout.readline, b''):
        chromosome, start, end, _ = stdout_line.decode().strip().split()
        yield [chromosome, int(start), int(end)]
    process.stdout.close()
    process.wait()


def _str_to_bed(bed_string):
    """Convert BED in string format into (str, int, int)."""
    chromosome, start, end = bed_string.strip().split("\t")
    start = int(start)
    end = int(end)
    return (chromosome, start, end)


def _get_fasta(transcriptome_dict, iterable_of_bed):
    """Extract subsequences in trancriptome_fn according to locis.

    (fasta file, list of lists) -> seqrecord
    """
    for bed in iterable_of_bed:
        chromosome, start, end = bed
        seq = transcriptome_dict[chromosome][start:end].seq
        identifier = "{0}:{1}-{2}".format(chromosome, start, end)
        description = identifier
        yield SeqRecord(id=identifier, seq=seq, description=description)


def _find_exons_pipeline(kmer, bloom_filter_fn, transcriptome_fn, max_fp_bases=5):
    """Find exons according to the Bloom filter -> BED
    Main pipeline:
    - Check every kmer,
    - merge if there is an overlap of k-1 bases
    - Throw away too short exons (k + max_fp_bases)
    - merge consecutive exons if they have an ovelap of 2*max_fp_bases
    """
    # Prepare the commands
    c_kmers = ["abyss-bloom", "kmers", "--kmer", str(kmer), "--verbose",
        "--bed", bloom_filter_fn, transcriptome_fn]
    c_merge1 = ["bedtools", "merge", "-d", str(-kmer + 2)]
    c_filter = ["awk", "$3 - $2 >= {min_length}".format(
        min_length=kmer + max_fp_bases
    )]
    c_merge2 = ["bedtools", "merge", "-d", str(-10 + 2)]
    # Run all of them streamlined
    p_kmers = Popen(c_kmers, stdout=PIPE)
    p_merge1 = Popen(c_merge1, stdin=p_kmers.stdout, stdout=PIPE)
    p_filter = Popen(c_filter, stdin=p_merge1.stdout, stdout=PIPE)
    p_merge2 = Popen(c_merge2, stdin=p_filter.stdout, stdout=PIPE)
    p_kmers.stdout.close()
    yield from map(_str_to_bed, _process_output(p_merge2))


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
    # Get predicted exons in bed format
    positive_exons_bed = _find_exons_pipeline(
        kmer, bloom_filter_fn, transcriptome_fn, max_fp_bases=5
    )
    # Bed -> fasta
    transcriptome_dict = SeqIO.index(filename=transcriptome_fn, format="fasta")
    positive_exons_fasta = _get_fasta(transcriptome_dict, positive_exons_bed)
    # Reduce
    exons = reduce_exons(positive_exons_fasta)  # Collapse identical exons into one
    # Write
    SeqIO.write(sequences=exons, handle=output_fasta, format="fasta")
