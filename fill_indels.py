#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(
    usage = 'python3 fill_indels.py -t treseq.tsv -i bloom.filter -k 27 ',
    description='Fix small indels within introns with help of abyss-sealer',
    epilog = 'Jorge Langa. Send issues and pull requests to github.com/jlanga/'
        'exon_finder'
)

parser.add_argument(
    '--input-treseq',
    '-t',
    type = str,
    required = True,
    help = 'Transcript Exon Sequence table from find_exons. Use - for stdin',
    dest = 'input_treseq',
    metavar = 'TSV'
)

parser.add_argument(
    '--input-bloom',
    '-i',
    type = str,
    required = True,
    help = 'Bloom filter with genomic sequences (from abyss-bloom)',
    dest = 'bloom_filter',
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

parser.add_argument(
    '--max-gap-size',
    '-G',
    type = int,
    required = False,
    help = 'The maximum size of the gap to be filled',
    dest = 'max_gap_size',
    default = 50,
    metavar = 'INT'
)

parser.add_argument(
    '--output-treseq',
    '-o',
    type = str,
    required = False,
    help = 'Path to output file. Default is stdout',
    dest = "output_treseq",
    metavar = "FILE"
)

args = parser.parse_args()



def treseq_to_fasta(treseq_file, fasta_file, ns = 100):
    """(str, str, int) -> new_file
    Open a treseq file and write it down into a fasta file. Exons belonging to 
    the same transcript are gapped in the same order provided. Space between
    two exons are filled with the number of Ns provided in ns.
    """
    from Bio import SeqIO
    from Bio import Seq
    transcript_seq = {}

    with open(treseq_file, "r") as treseq_in, open(fasta_file, "w") as fasta_out:

        # Dump all data into a dict where key is transcript id and key is the gapped sequence
        # tr121312 : CACACTAGCTAGCTANNNNNNNNNNNNNNNNNNNNNACATCGATCGTACGTGACTGANNNNNNNNNNCACAC
        for line in treseq_in:
            (transcript_id, exon_id, sequence) = line.split()

            if transcript_id not in transcript_seq:
                transcript_seq[transcript_id] = sequence
            else:
                transcript_seq[transcript_id] += 'N' * ns + sequence

        # Write all data into fasta
        for transcript_id, sequence in transcript_seq.items():
            fasta_out.write(
                ">" + transcript_id + "\n" +
                sequence + "\n"
            )





def fasta_to_treseq(fasta_file, treseq_file):
    """(str, str) -> file
    Transform a fasta file with records with Ns into a treseq file, i.e., given a fasta record
    of the form:
    >seq1
    AGTCAGCTAGCANNNNNNCAGTCAGCTGACNNNNNNCGATCGATCGA

    the result is
    seq1    seq1_e1 AGTCAGCTAGCA
    seq1    seq1_e2 CAGTCAGCTGAC
    seq1    seq1_e3 CGATCGATCGA
    """
    from Bio import SeqIO

    with open(fasta_file, 'r') as fasta_in, open(treseq_file, 'w') as treseq_out:
        records = SeqIO.parse(format= "fasta", handle = fasta_in)
        for record in records:

            # Extract id and seq
            identifier = record.id
            sequence = str(record.seq)

            # Split transcripts into exons based on the number of Ns
            exons = tuple(x for x in str(sequence).split('N') if x != '')

            # Write to file
            for exon_number, exon_sequence in enumerate(exons):

                exon_name = identifier + "_e" + str(exon_number + 1)

                treseq_out.write("\t".join([
                    identifier, exon_name, exon_sequence
                ]) + "\n")



def run_sealer(input_fasta_fn, kmer, bloom_filter_fn, output_prefix, max_gap_size=50):
    """(str, int, str, str, int) -> file
    Run abyss-sealer with:
        - gapped fasta file from input_fasta_fn
        - kmer size kmer
        - bloom filter already constructed from bloom_filter_fn,
        - output_prefix as prefix for scaffold, merge and log file
        - a max gap size of max_gap_size
    """
    from subprocess import Popen, PIPE
    from sys import stderr
    
    # Compose sealer command
    sealer_command = [  # Run abyss-sealer kmers
        "abyss-sealer",
            "--input-scaffold", input_fasta_fn,
            "--flank-length", str(kmer),
            "--max-gap-length", str(max_gap_size),
            "--kmer", str(kmer),
            "--input-bloom", bloom_filter,
            "--fix-errors",
            "--mask",
            "--search-mem", "1G",
            "--output-prefix", output_prefix,
            "--verbose"
    ]
    print(sealer_command)
    # Run the pipeline
    sealer_process = Popen(sealer_command, stderr= PIPE)
   

    # Throw away sdterr
    for line in sealer_process.stderr.readlines():
        stderr.write("[abyss-sealer]\t" + line.decode())





if __name__ == "__main__":
    
    
    # Read args
    args = vars(parser.parse_args())
    input_treseq = args["input_treseq"]
    bloom_filter = args["bloom_filter"]
    kmer = args["kmer"]
    max_gap_size = args["max_gap_size"]
    output_treseq = args["output_treseq"]

    # Prepare temporary files
    from tempfile import NamedTemporaryFile
    prefix = NamedTemporaryFile()
    fasta_fn = prefix.name + ".fa"
    scaffold_fn = prefix.name + "_scaffold.fa"
    merged_fn = prefix.name + "_merged.fa"
    log_fn = prefix.name + "_log.txt"

    # treseq_to_fasta
    treseq_to_fasta(input_treseq, fasta_fn, 100)
    # run sealer over fasta
    run_sealer(
        input_fasta_fn= fasta_fn,
        kmer= kmer,
        bloom_filter_fn= bloom_filter,
        output_prefix= prefix.name,
        max_gap_size= max_gap_size
    )

    # fasta to treseq
    fasta_to_treseq(scaffold_fn, output_treseq)

    # Delete temp files
    import os
    prefix.close()
    os.remove(fasta_fn)
    os.remove(scaffold_fn)
    os.remove(merged_fn)
    os.remove(log_fn)