#!/usr/bin/env python3

# Import everything
from subprocess import Popen, PIPE
from Bio.SeqRecord import SeqRecord


def _process_output(process):
    """Get lines in bed format from the output of a Popen."""
    for stdout_line in iter(process.stdout.readline, b''):
        chromosome, start, end = stdout_line.decode().strip().split()
        yield (chromosome, int(start), int(end))
    process.stdout.close()
    process.wait()


def _get_fasta(transcriptome_dict, iterable_of_bed):
    """Extract subsequences in trancriptome_fn according to locis.

    (fasta file, list of lists) -> seqrecord
    """
    for bed in iterable_of_bed:
        chromosome, start, end = bed
        if chromosome in transcriptome_dict:
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
    yield from _process_output(p_merge2)
