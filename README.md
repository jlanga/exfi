[![Build Status](https://travis-ci.org/jlanga/exfi.svg?branch=master)](https://travis-ci.org/jlanga/exfi)

# exfi
Get exons from a transcriptome and raw genomic reads using abyss-bloom and bedtools

## Requirements
```
abyss>=2.0.0
bedtools (tested on 2.0)
python3
biopython
```
We recomend installing these packages with `conda` and `bioconda`

## How to install

Copy this repo and install it with `pip`:

```sh
git clone https://github.com/jlanga/exfi.git
pip install --user exfi
```

## Required data

- A transcriptome in fasta format (take it from Ensembl for example, or the result of a _de novo_ transcriptome assembler like trinity, trans-abyss or oases)

- A set of genomic reads in fastq format, paired end or not. `.gz` files are allowed.

## Running the pipeline

1. Make a Bloom filter of the genomic reads with `abyss-bloom build`. 
- `genome.fa.gz` is the set of genomic reads and
- `genome_k27_m500M_l1.bloom` is the resulting Bloom filter, made of kmers of length 27, a size of 500 Mb and the number of times of a kmer must be in the reads is 1 (levels).

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

2. Run `find_exons` to get putative exons in the transcriptome.
- `data/transcript.fa` is the input transcriptome,
- `genome_k27_m500M_l1.bloom` is the Bloom filter generated above
- kmer length has to be the same
- `test_exons_raw.fa` is the first set of putative kmers.
```sh
find_exons \
    --input-fasta data/transcript.fa \
    --input-bloom genome_k27_m500M_l1.bloom \
    --kmer 27 \
    --output-fasta test_exons_raw.fa
```


3.  Filter the exons

The raw set of exons is filled with the false positives that the Bloom filter introduces, made of 
- exons of length ~k
- false bases at the ends of the true exons
Therefore it is necessary to clean it. We provide two methods:

3.1 By extensibility
Try if the kmers at the beginning and the end of the exon are extendable, i.e., exist at least a 1-base extension of those two extremes:

```sh
filter_exons_by_extensibility \
    --input-fasta test_exons_raw.fa \
    --input-bloom genome_k27_m500M_l1.bloom \
    --kmer 27 \
    --output-fasta test_exons_filtered_by_extensiblity.fa
```

3.2 By length
Cut from each of the exons one or two bases, and drop the exons that have a length smaller than a certain length:
```sh
filter_exons_by_length \
    --input-fasta test_exons_raw.fa \
    --minimum-length 27 \
    --trim-left 1 \
    --trim-right 1 \
    --output-fasta test_exons_filtered_by_length.fa
```

## Results

The output is a fasta file in which the sequences are the different exons found. The header contains a unique identifier and the coordinates where the exons can be found in the transcriptome.

```
>EXONNUMBER TRANSCRIPT_ID1:start-end TRANSCRIPT_ID2:start-end
```

```
>EXON00000000022 ENSDART00000163017.1:419-657 ENSDART00000172567.1:406-644
CAGGAGGAAGATGGCCTTACGATGCTGCAGTTGGAGAAGGACCTCCGCACGCAGGTCAAG
TTGATGCTTAAAGAAAAGAGCTCTCGCCAGACTGAGCTGAAGTCTCTGATTCAGCAGGAT
CAAGATCTGTGTGATGTCCTGTGTGAAGACCTGTTCCCCATCCATCCCGAACGCGTGCCT
TCACAGCAGCAGCTGCAGAACTACAGGCAGCACATCAACACACGCAACCAGGAGAAG
```

## Attributions

Written by 

## Bibliography

- [abyss](https://github.com/bcgsc/abyss/)

- [bedtools](https://bedtools.readthedocs.io/)

- [biopython](biopython.org)
