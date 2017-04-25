[![Build Status](https://travis-ci.org/jlanga/exon_finder.svg?branch=master)](https://travis-ci.org/jlanga/exon_finder)

# exon_finder
Get exons from a transcriptome and raw genomic reads using abyss-bloom and bedtools

## Requirements
```
abyss>=2.0.0
bedtools (tested on 2.0)
python3
biopython
```
We recomend installing these packages with `conda` and `bioconda`

## Running the pipeline

1. Make a Bloom filter of the genomic reads with `abyss-bloom build`:

```sh
abyss-bloom build \
    --kmer=27 \
    --verbose \
    --bloom-size=500M \
    --threads=4 \
    --levels=1 \
    genome_k27_m500M_l1.bloom \
    test/genome.fa.gz
```

2. Run `find_exons` to get putative exons in the transcriptome:
```sh
find_exons
    --input-fasta data/transcript.fa
    --input-bloom genome_k27_m500M_l1.bloom
    --kmer 27
    --output-fasta test_exons_raw.fa
```


3.  Filter the exons by extensibility:
```sh
filter_exons_by_extensibility
    --input-fasta test_exons_raw.fa
    --input-bloom genome_k27_m500M_l1.bloom
    --kmer 27
    --output-fasta test_exons_filtered_by_extensiblity.fa
```


## Results

The output is a fasta file in which the sequences are the different exons found. The header info is as follows:

```
>EXONNUMBER TRANSCRIPT_ID1:start-end;rank TRANSCRIPT_ID2:start-end;rank
```

```
>EXON00000000022 ENSDART00000163017.1:419-657;4 ENSDART00000172567.1:406-644;4
CAGGAGGAAGATGGCCTTACGATGCTGCAGTTGGAGAAGGACCTCCGCACGCAGGTCAAG
TTGATGCTTAAAGAAAAGAGCTCTCGCCAGACTGAGCTGAAGTCTCTGATTCAGCAGGAT
CAAGATCTGTGTGATGTCCTGTGTGAAGACCTGTTCCCCATCCATCCCGAACGCGTGCCT
TCACAGCAGCAGCTGCAGAACTACAGGCAGCACATCAACACACGCAACCAGGAGAAG
```
