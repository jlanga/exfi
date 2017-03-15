#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(
    usage = 'python3 exon_finder.py -t transcriptome.fa -i bloom.filter -k 27 '
            '-o gene.fa',
    description='Process a transcriptome to extract exonic regions',
    epilog = 'Jorge Langa. Send issues and pull requests to github.com/jlanga/'
        'exon_finder'
)

parser.add_argument(
    '--input-transcriptome',
    '-t',
    type = str,
    required = True,
    help = 'Transcriptome from which to extract the exons',
    dest = 'transcriptome',
    metavar = 'FILE_FASTA'
)

parser.add_argument(
    '--input-bloom',
    '-i',
    type = str,
    required = True,
    help = 'Bloom filter with genomic sequences (from abyss-bloom)',
    dest = 'bloom',
    metavar = 'BLOOM'
)

parser.add_argument(
    '--kmer',
    '-k',
    type = int,
    required = True,
    help = 'The size of the k-mer',
    dest = 'kmer',
    metavar = 'KMER'
)

args = parser.parse_args()


def run_pipeline(transcriptome_fn, kmer, bloom_filter_fn):
    
    from subprocess import Popen, PIPE
    from sys import stdin, stdout, stderr
    
    command1 = [  # Run abyss-bloom kmers
    "abyss-bloom", "kmers",
        "--kmer", str(kmer),
        "--verbose",
        "--bed",
        bloom_filter_fn,
        transcriptome_fn
    ]

    command2 = [  # Merge overlapping kmers
    "bedtools", "merge",
        "-d",
        str(- kmer + 2)
    ]

    command3 = [  # Filter lonely kmers (most likely False positives)
        'awk', 
            " ".join(['"{if($3 - $2 > ',  str(kmer), ') print}"'])
    ]

    command4 = [  # Get transcriptid:coordinates TAB sequence
        "bedtools", "getfasta",
            "-fi", transcriptome_fn,
            "-bed", "-",
            "-tab"
    ]
    
    # Run the pipeline
    p1 = Popen(command1, stdout= PIPE, stderr= PIPE)
    p2 = Popen(command2, stdin=p1.stdout, stdout= PIPE, stderr= PIPE)
    p3 = Popen(command3, stdin=p2.stdout, stdout= PIPE, stderr= PIPE)
    p4 = Popen(command4, stdin=p3.stdout, stdout= PIPE, stderr= PIPE)
    
    # Throw away 
    for line in p1.stderr.readlines():
        stderr.write(line.decode())
    for line in p2.stderr.readlines():
        stderr.write(line.decode())
    for line in p3.stderr.readlines():
        stderr.write(line.decode())
    for line in p4.stderr.readlines():
        stderr.write(line.decode())
    
    # Reformat output
    ## Dict with the exons per transcript seen
    transcript_to_exon_number = {}
    
    for line in p4.stdout.readlines():
        line = line.decode()
        [transcript_id, sequence] = line.split("\t")
        transcript_id = transcript_id.rsplit(":")[0]
        if transcript_id in transcript_to_exon_number:
            transcript_to_exon_number[transcript_id] += 1
        else:
            transcript_to_exon_number[transcript_id] = 1
        exon_number = transcript_to_exon_number[transcript_id]
        exon_id = transcript_id + "_exon"+ str(exon_number)
        stdout.write(
           "\t".join([transcript_id, exon_id, sequence])
        )

if __name__ == "__main__":
    
    import sys
    
    args = vars(parser.parse_args())
    transcriptome = args["transcriptome"]
    bloom_filter = args["bloom"]
    kmer = args["kmer"]

    run_pipeline(transcriptome, kmer, bloom_filter)