#!/usr/bin/env python3

# Import everything
from subprocess import Popen, PIPE
from sys import stdin, stdout, stderr
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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

# def _bedtools_merge_command(kmer):
#     """Merge overlapping kmers"""
#     return ["bedtools", "merge",
#         "-d", str(- kmer + 2)
#     ]
#
# def _bedtools_getfasta_command(transcriptome_fn):
#     """Get transcriptid:coordinates TAB sequence"""
#     return [ "bedtools", "getfasta",
#         "-fi", transcriptome_fn,
#         "-bed", "-",
#         "-tab"
#     ]

def _process_output(process):
    """
    Auxiliar function to process the output as it comes

    NEEDS IMPROVEMENT
    """
    raw = process.communicate()[0]
    decoded = raw.decode()
    #del raw
    splitted = decoded.split("\n")
    #del decoded
    for line in splitted:
        yield line

def _abyss_kmers_to_bed(iterable_of_str):
    """
    Convert BED in string format into str, int, int
    """
    for record in iterable_of_str:
        if record:
            chromosome, start, end, _ = record.strip().split()
            start = int(start)
            end = int(end)
            yield [chromosome, start, end]

def _merge_bed(bed_records):
    """
    Merge overlapping by all but one base records
    """
    parsed_records = _abyss_kmers_to_bed(bed_records)
    old = [None, None, None]  # Loci, start, end
    for new in parsed_records:
        if not old[0]: # First record
            old = new
            continue
        if old[0] != new[0]: # Change of loci
            yield old
            old = new
            continue
        if old[0] == new[0] and new[2] != old[2] +1: # Same record
            yield old
            old = new
            continue
        old[2] = new[2]  # If not different contig or too big jump, update
    yield old # last step

def _get_fasta(transcriptome_fn, locis):
    transcriptome_dict = SeqIO.to_dict(SeqIO.parse(transcriptome_fn, "fasta"))
    for loci in locis:
        if not loci:
            continue
        chromosome, start, end = loci #.strip().split("\t")
        seq = transcriptome_dict[chromosome][int(start):int(end)].seq
        identifier = "{0}:{1}-{2}".format(chromosome, start, end)
        description = identifier
        #print(record)
        yield SeqRecord(
            id= identifier,
            seq= seq,
            description = description
        )


def _find_exons_pipeline(kmer, bloom_filter_fn, transcriptome_fn):
    # Prepare the commands
    abyss_bloom_kmers = _abyss_bloom_kmers_command(
        kmer, bloom_filter_fn, transcriptome_fn
    )
    #bedtools_merge = _bedtools_merge_command(kmer)
    #bedtools_getfasta = _bedtools_getfasta_command(transcriptome_fn)

    p1 = Popen(abyss_bloom_kmers, stdout= PIPE)
    #p2 = Popen(bedtools_merge, stdin= p1.stdout, stdout= PIPE)

    # Manage all processes properly
    #p1.stdout.close()

    abyss_kmers_output = _process_output(p1)
    merged = _merge_bed(abyss_kmers_output)
    records = _get_fasta(transcriptome_fn, merged)

    for line in records:
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
    #seqrecords = tab_to_seqrecord(pipeline_output)
    exons = reduce_exons(pipeline_output)  # Collapse identical exons into one
    SeqIO.write(
        sequences= exons,
        handle= output_fasta,
        format= "fasta"
    )
