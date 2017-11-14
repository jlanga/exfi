[![Build Status](https://travis-ci.org/jlanga/exfi.svg?branch=master)](https://travis-ci.org/jlanga/exfi)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/e3200fef4f7549d78c2cf85364f1c602)](https://www.codacy.com/app/jorge.langa.arranz/exfi?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=jlanga/exfi&amp;utm_campaign=Badge_Grade)
[![Codacy Badge](https://api.codacy.com/project/badge/Coverage/e3200fef4f7549d78c2cf85364f1c602)](https://www.codacy.com/app/jorge.langa.arranz/exfi?utm_source=github.com&utm_medium=referral&utm_content=jlanga/exfi&utm_campaign=Badge_Coverage)
# exfi
Get exons from a transcriptome and raw genomic reads using abyss-bloom and bedtools

## Requirements
```
abyss>=2.0.0
bedtools (tested on 2.0)
python3
biopython
networkx
pandas
biobloomtools
```
We recomend installing these packages with `conda` and `bioconda`

## How to install

Copy this repo and install it with `pip`:

```sh
git clone https://github.com/jlanga/exfi.git
pip install --user exfi
```

To install other dependencies, follow the instructions from the travis files:

1. Install packages with `apt`:

```sh
sudo apt install build-essential git curl libboost-dev gcc autoconf bzip2 zlib1g libsparsehash-dev
```

2. Install conda, then

```sh
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda install --yes abyss biopython bedtools networkx pandas pip
```

3. Install `biobloomtools`

There is an easy way, with `brew` (not the latest release):

```sh
brew install biobloomtools
```

Or the manual way:

```sh
# Install SDSL-Lite
git clone https://github.com/simongog/sdsl-lite.git
pushd sdsl-lite/ && \
sudo ./install.sh /usr/local/ && \
popd

# Install biobloomtools
git clone https://github.com/bcgsc/biobloom.git
pushd biobloom/ && \
git submodule update --init && \
./autogen.sh && \
./configure --prefix=/usr/local/ && \
make -j 4 && \
sudo make install && \
popd
```



## Required data

- A transcriptome in fasta format (take it from Ensembl for example, or the result of a _de novo_ transcriptome assembler like trinity, trans-abyss or oases)

- A set of genomic reads in fastq format, paired end or not. `.gz` files are allowed.

## Running the pipeline

1. Make a baited Bloom filter of the genomic reads with `build_baited_bloom_filter`:
- `genome.fa.gz` is the set of genomic reads and
- `genome_k27_m500M_l1.bloom` is the resulting Bloom filter, made of kmers of length 27, a size of 500 Mb and the number of times of a kmer must be in the reads is 1 (levels).

```sh
build_baited_bloom_filter \
    --input-fasta data/transcript.fa \
    --kmer 27 \
    --bloom-size 500M \
    --levels 1 \
    --threads 4 \
    --output-bloom results/genome_k27_m500M_l1.bloom \
    genome.fa.gz
```

2. Run `build_splice_graph` to get putative exons in the transcriptome.
- `data/transcript.fa` is the input transcriptome,
- `genome_k27_m500M_l1.bloom` is the Bloom filter generated above
- kmer length has to be the same
- `test.gfa` is the resulting splice graph in [GFA1 format](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md).

```sh
build_splicegraph \
    --input-fasta data/transcript.fa \
    --input-bloom results/genome_k27_m500M_l1.bloom \
    --kmer 27 \
    --max-fp-bases 5 \
    --output-gfa test.gfa
```

This splice graph can be visualized with [Bandage](https://rrwick.github.io/Bandage/)

Example:

```
H       VN:Z:1.0
S       EXON00000000001 GTAAGCCGCGGCGGTGTGTGTGTGTGTGTGTGTTCTCCGTCATCTGTGTTCTGCTGAATGATGAGGACAGACGTGTTTCTCCAGCGGAGGAAGCGTAGAGATGTTCTGCTCTCCATCATCGCTCTTCTTCTGCTCATCTTCGCCATCGTTCATCTCGTCTTCTGCGCTGGACTGAGTTTCCAGGGTTCGAGTTCTGCTCGCGTCCGCCGAGACCTC        LN:i:216
S       EXON00000000002 GAGAATGCGAGTGAGTGTGTGCAGCCACAGTCGTCTGAGTTTCCTGAAGGATTCTTCACGGTGCAGGAGAGGAAAGATGGAGGAATCCTGATTTACTTCATGATCATCTTCTACATGCTGCTGTCCGTCTCCATCGTGTGTGATGAATATTTTCTGCCATCTCTGGAGGTCATCAGCGAGCG  LN:i:182
S       EXON00000000003 GTCTTGGTCTCTCGCAGGATGTTGCTGGAGCCACGTTTATGGCTGCGGGGAGTTCGGCTCCAGAGCTCGTCACTGCATTTCTGGG   LN:i:85
S       EXON00000000004 GGTGTGTTTGTGACGAAGGGCGACATCGGCGTCAGCACCATCATGGGTTCTGCTGTCTATAACCTGCTGTGCATCTGTGCAGCGTGCGGCCTGCTGTCCTCTGCAG      LN:i:106
S       EXON00000000005 GTTGGTCGTCTGAGCTGCTGGCCGTTGTTCAGAGATTGTGTTGCGTACTCCATCAGTGTCGCCGCCGTCATCGCCATCATCTCAGATAACAGAGTTTACTGG  LN:i:102
S       EXON00000000006 GGTATGATGGCGCGTGTCTCCTGCTGGTGTACGGTGTGTATGTAGCTGTACTGTGTTTCGATCTGAAGATCAGCGAGTACGTGATGCAGCGCTTCAGTCCATGCTGCTGGTGTCTGAAACCTCGCGATCGTGACTCAGGCGAGCAGCAGCCTCTAGTGGGCTGGAGTGACGACAGCAGCCTGCGGGTCCAGCGCCGTTCCAGAAATGACAGCGGAATATTCCAGGATGATTCTGGATATTCACATCTATCGCTCAGCCTGCACGGACTCAACGAAATCAGCGAC    LN:i:284
S       EXON00000000007 GAGCACAAGAGTGTGTTCTCCATGCCGGATCACGATCTGAAGCGAATCCTGTGGGTTTTGTCTCTTCCGGTCAGCACTCTGCTGTTTGTGAGCGTTCCCGACTGCAGGAGACCCTTCTGGAAGAACTTCTACATGCTGACCTTCCTGATGTCCGCCGTCTGGATTTCTGCATTCACTTATGTGCTGGTCTGGATGGTCACAATCGTGG        LN:i:208
S       EXON00000000008 GGGGAGACTCTGGGAATCCCGGACACAGTGATGGGAATGACTCTTCTGGCTGCAGGAACCAGTATCCCCGACACCGTGGCCAGTGTGATGGTGGCACGAGAAGGTAA     LN:i:107
S       EXON00000000009 AGGTAAATCTGATATGGCCATGTCCAACATCGTGGGCTCTAACGTGTTCGATATGCTGTGTCTGGGCCTGCCGTGGTTCATCCAGACGGTGTTTGTTGACGTGGGCTCCCCGGTGGATGTCAACAGCTCGGGGCTGGTCTTCATGTCCTGCACGCTGCTGCTCTCCATCATCTTCCTCTTCCTCGCCGTGCACATCAACGGCTGGAAGCTGGACTGGAAGCTGGGTCTGGTGTGTTTGGCGTGTTACATTCTGTTCGCAACACTCTCCATCCTGTACGAGCTCGGCATCATCGGGAACAATCCCATACGCTCCTGCAGCGACTGAACACTGCTCTACAGCGCCCCCTTATGGACAACACAAGGACGTGACTCTTTATAACCCTCTAAAGTGCACAGGTTCATTACTGAATACAAGAAAATAGAACTGCGAGACGTCAACTCAAAATACAAGAGAAGTCAAAGTGCGAGATGTAAAAAATATATGCACATAAATGAGGATAAACTTTTTATTTAATAAGACAAAACTGCATAAAGTCTGATGTGAACACTGCTCAACAGCGCCCTCTCATGGACAACACATGGATCTGACTCTTATTAACCCTCCAGAGTGCAAATACACTAACACAACGTAATATAACCAAGTTAAAATGGCAAGATGTGAACTCAAAATACAAGAAAGCAGTCAAGATGCCCGACATAACAAATGTGCATTAAAATGTAAGCCC   LN:i:725
L       EXON00000000001 +       EXON00000000002 +       0M
L       EXON00000000002 +       EXON00000000003 +       1M
L       EXON00000000003 +       EXON00000000004 +       2M
L       EXON00000000004 +       EXON00000000005 +       1M
L       EXON00000000005 +       EXON00000000006 +       2M
L       EXON00000000006 +       EXON00000000007 +       0M
L       EXON00000000007 +       EXON00000000008 +       1M
L       EXON00000000008 +       EXON00000000009 +       6M
C       ENSDART00000033574.5    +       EXON00000000001 +       0       216M
C       ENSDART00000033574.5    +       EXON00000000002 +       216     182M
C       ENSDART00000033574.5    +       EXON00000000003 +       397     85M
C       ENSDART00000033574.5    +       EXON00000000004 +       480     106M
C       ENSDART00000033574.5    +       EXON00000000005 +       585     102M
C       ENSDART00000033574.5    +       EXON00000000006 +       685     284M
C       ENSDART00000033574.5    +       EXON00000000007 +       969     208M
C       ENSDART00000033574.5    +       EXON00000000008 +       1176    107M
C       ENSDART00000033574.5    +       EXON00000000009 +       1277    725M
P       ENSDART00000033574.5    EXON00000000001+,EXON00000000002+,EXON00000000003+,EXON00000000004+,EXON00000000005+,EXON00000000006+,EXON00000000007+,EXON00000000008+,EXON00000000009+

```

3.  Get exonic sequences

To extract meaningful information from the GFA file, we provide two scripts:

- `gfa_to_exons`: which returns the predicted exons in FASTA format. For each record, each sequence comes with a unique identifier (`EXON[0-9]+`), a description indicating the coordinates of this exon with respect to the different transcripts (`TR1:0-200 TR2:105-305`), and the sequence of nucleotides. It is possible to hard and soft mask nucleotides that may not be correct. Example (soft masked):

```
>EXON00000000003 ENSDART00000033574.5:397-482
gTCTTGGTCTCTCGCAGGATGTTGCTGGAGCCACGTTTATGGCTGCGGGGAGTTCGGCTC
CAGAGCTCGTCACTGCATTTCTGgg
>EXON00000000004 ENSDART00000033574.5:480-586
ggTGTGTTTGTGACGAAGGGCGACATCGGCGTCAGCACCATCATGGGTTCTGCTGTCTAT
AACCTGCTGTGCATCTGTGCAGCGTGCGGCCTGCTGTCCTCTGCAg
```

- `gfa_to_gapped_transcript`: which returns the transcript with interleaved `N`s where it is predicted to be an intron. Example (hard masked):
```
>ENSDART00000033574.5 EXON00000000001,EXON00000000002,EXON00000000003,EXON00000000004,EXON00000000005,EXON00000000006,EXON00000000007,EXON00000000008,EXON00000000009
GTAAGCCGCGGCGGTGTGTGTGTGTGTGTGTGTTCTCCGTCATCTGTGTTCTGCTGAATG
ATGAGGACAGACGTGTTTCTCCAGCGGAGGAAGCGTAGAGATGTTCTGCTCTCCATCATC
GCTCTTCTTCTGCTCATCTTCGCCATCGTTCATCTCGTCTTCTGCGCTGGACTGAGTTTC
CAGGGTTCGAGTTCTGCTCGCGTCCGCCGAGACCTCNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNGAGAATGCGAGTGAGTGTGTGCAGCCACAGTCGTCTGAGTTTCC
TGAAGGATTCTTCACGGTGCAGGAGAGGAAAGATGGAGGAATCCTGATTTACTTCATGAT
CATCTTCTACATGCTGCTGTCCGTCTCCATCGTGTGTGATGAATATTTTCTGCCATCTCT
GGAGGTCATCAGCGAGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNT
CTTGGTCTCTCGCAGGATGTTGCTGGAGCCACGTTTATGGCTGCGGGGAGTTCGGCTCCA
GAGCTCGTCACTGCATTTCTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNTGTGTTTGTGACGAAGGGCGACATCGGCGTCAGCACCATCATGGGTTCTGCTGTC
TATAACCTGCTGTGCATCTGTGCAGCGTGCGGCCTGCTGTCCTCTGCANNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTGGTCGTCTGAGCTGCTGGCCGTTGTTCA
GAGATTGTGTTGCGTACTCCATCAGTGTCGCCGCCGTCATCGCCATCATCTCAGATAACA
GAGTTTACTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTATGATG
GCGCGTGTCTCCTGCTGGTGTACGGTGTGTATGTAGCTGTACTGTGTTTCGATCTGAAGA
TCAGCGAGTACGTGATGCAGCGCTTCAGTCCATGCTGCTGGTGTCTGAAACCTCGCGATC
GTGACTCAGGCGAGCAGCAGCCTCTAGTGGGCTGGAGTGACGACAGCAGCCTGCGGGTCC
AGCGCCGTTCCAGAAATGACAGCGGAATATTCCAGGATGATTCTGGATATTCACATCTAT
CGCTCAGCCTGCACGGACTCAACGAAATCAGCGACNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNGAGCACAAGAGTGTGTTCTCCATGCCGGATCACGATCTGAAGCGA
ATCCTGTGGGTTTTGTCTCTTCCGGTCAGCACTCTGCTGTTTGTGAGCGTTCCCGACTGC
AGGAGACCCTTCTGGAAGAACTTCTACATGCTGACCTTCCTGATGTCCGCCGTCTGGATT
TCTGCATTCACTTATGTGCTGGTCTGGATGGTCACAATCGTGNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNGGGAGACTCTGGGAATCCCGGACACAGTGATGGGAA
TGACTCTTCTGGCTGCAGGAACCAGTATCCCCGACACCGTGGCCAGTGTGATGGTGGCAC
GAGANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNATCT
GATATGGCCATGTCCAACATCGTGGGCTCTAACGTGTTCGATATGCTGTGTCTGGGCCTG
CCGTGGTTCATCCAGACGGTGTTTGTTGACGTGGGCTCCCCGGTGGATGTCAACAGCTCG
GGGCTGGTCTTCATGTCCTGCACGCTGCTGCTCTCCATCATCTTCCTCTTCCTCGCCGTG
CACATCAACGGCTGGAAGCTGGACTGGAAGCTGGGTCTGGTGTGTTTGGCGTGTTACATT
CTGTTCGCAACACTCTCCATCCTGTACGAGCTCGGCATCATCGGGAACAATCCCATACGC
TCCTGCAGCGACTGAACACTGCTCTACAGCGCCCCCTTATGGACAACACAAGGACGTGAC
TCTTTATAACCCTCTAAAGTGCACAGGTTCATTACTGAATACAAGAAAATAGAACTGCGA
GACGTCAACTCAAAATACAAGAGAAGTCAAAGTGCGAGATGTAAAAAATATATGCACATA
AATGAGGATAAACTTTTTATTTAATAAGACAAAACTGCATAAAGTCTGATGTGAACACTG
CTCAACAGCGCCCTCTCATGGACAACACATGGATCTGACTCTTATTAACCCTCCAGAGTG
CAAATACACTAACACAACGTAATATAACCAAGTTAAAATGGCAAGATGTGAACTCAAAAT
ACAAGAAAGCAGTCAAGATGCCCGACATAACAAATGTGCATTAAAATGTAAGCCC
```


## Attributions

Written by Jorge Langa

## Bibliography

- [abyss](https://github.com/bcgsc/abyss/)

- [bandage](https://rrwick.github.io/Bandage/)

- [bedtools](https://bedtools.readthedocs.io/)

- [biobloomtools](https://github.com/bcgsc/biobloom)

- [biopython](http://biopython.org/)

- [networkx](https://networkx.github.io/)

- [pandas](http://pandas.pydata.org/)

- [funniest](https://pypi.python.org/pypi/funniest/0.1)

- [sphinx](www.sphinx-doc.org)

- [unittest](https://docs.python.org/3/library/unittest.html)
