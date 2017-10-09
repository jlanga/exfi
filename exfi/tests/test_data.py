#!/usr/bin/env python3

from exfi.io import \
    _clean_index

from Bio import SeqIO

path_simple = {"ENSDART00000161035.1": ["EXON00000000001"]}

path_different = {
    "ENSDART00000161035.1":
        ["EXON00000000001", "EXON00000000002", "EXON00000000003"]
    ,
    "ENSDART00000165342.1":
        ["EXON00000000004", "EXON00000000005", "EXON00000000006",
        "EXON00000000007", "EXON00000000008", "EXON00000000009",
        "EXON00000000010", "EXON00000000011", "EXON00000000012",
        "EXON00000000013", "EXON00000000014", "EXON00000000015"]
}

index_simple = _clean_index(SeqIO.index(
    filename="exfi/tests/files/build_splicegraph/single.fa",
    format="fasta"
))

index_different = _clean_index(SeqIO.index(
    filename="exfi/tests/files/build_splicegraph/different_transcripts.fa",
    format="fasta"
))

transcriptome_simple = _clean_index(SeqIO.index(
    filename="exfi/tests/files/build_splicegraph/transcriptome_simple.fa",
    format="fasta"
))

transcriptome_different = _clean_index(SeqIO.index(
    filename="exfi/tests/files/build_splicegraph/transcriptome_different.fa",
    format="fasta"
))


empty_gfa = "exfi/tests/files/io/empty.gfa"
single_gfa = "exfi/tests/files/io/single.gfa"
different_gfa = "exfi/tests/files/io/different.gfa"

empty_exons = "exfi/tests/files/io/empty_exons.fa"
single_exons = "exfi/tests/files/io/single_exons.fa"
different_exons = "exfi/tests/files/io/different_exons.fa"
different_exons_masked = "exfi/tests/files/io/different_exons_masked.fa"

empty_gapped = "exfi/tests/files/io/empty_gapped.fa"
single_gapped = "exfi/tests/files/io/single_gapped.fa"
different_gapped = "exfi/tests/files/io/different_gapped.fa"
different_gapped_masked = "exfi/tests/files/io/different_gapped_masked.fa"
