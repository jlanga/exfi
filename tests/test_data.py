#!/usr/bin/env python3

from exfi.io import \
    index_fasta

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

index_simple = index_fasta(
    filename="tests/files/build_splicegraph/single.fa",
)

index_different = index_fasta(
    filename="tests/files/build_splicegraph/different_transcripts.fa",
)

transcriptome_simple = index_fasta(
    filename="tests/files/build_splicegraph/transcriptome_simple.fa",
)

transcriptome_different = index_fasta(
    filename="tests/files/build_splicegraph/transcriptome_different.fa",
)


empty_gfa = "tests/files/io/empty.gfa"
single_gfa = "tests/files/io/single.gfa"
different_gfa = "tests/files/io/different.gfa"

empty_exons = "tests/files/io/empty_exons.fa"
single_exons = "tests/files/io/single_exons.fa"
different_exons = "tests/files/io/different_exons.fa"
different_exons_soft = "tests/files/io/different_exons_soft.fa"
different_exons_hard = "tests/files/io/different_exons_hard.fa"

empty_gapped = "tests/files/io/empty_gapped.fa"
single_gapped = "tests/files/io/single_gapped.fa"
different_gapped = "tests/files/io/different_gapped.fa"
different_gapped_soft = "tests/files/io/different_gapped_soft.fa"
different_gapped_hard = "tests/files/io/different_gapped_hard.fa"
