#!/usr/bin/env bash
set -euxo pipefail


mkdir -p results/

if [ ! -f results/test_k27_m500M_l1.bloom ]; then
    build_baited_bloom_filter \
        --input-fasta data/transcript.fa \
        --kmer-size 27 \
        --bloom-size 500M \
        --levels 1 \
        --threads 4 \
        --output-bloom results/test_k27_m500M_l1.bloom \
        --verbose \
        data/genome.fa.gz
fi

build_splice_graph \
    --input-fasta data/transcript.fa \
    --input-bloom results/test_k27_m500M_l1.bloom \
    --kmer-size 27 \
    --max-fp-bases 5 \
    --output-gfa results/test.gfa \
    --verbose \
    --threads 4

build_splice_graph \
    --input-fasta data/transcript.fa \
    --input-bloom results/test_k27_m500M_l1.bloom \
    --kmer-size 27 \
    --max-fp-bases 5 \
    --output-gfa results/test_correct.gfa \
    --verbose \
    --threads 4 \
    --correct

build_splice_graph \
    --input-fasta data/transcript.fa \
    --input-bloom results/test_k27_m500M_l1.bloom \
    --kmer-size 27 \
    --max-fp-bases 5 \
    --output-gfa results/test_correct_polish.gfa \
    --verbose \
    --threads 4 \
    --correct \
    --polish

gfa1_to_fasta \
    --input-gfa results/test.gfa \
    --output-fasta results/test_exons.fa \
    --soft-mask-overlaps \
    --verbose

gfa1_to_fasta \
    --input-gfa results/test.gfa \
    --output-fasta results/test_gapped.fa \
    --hard-mask-overlaps \
    --gapped-transcript \
    --verbose
