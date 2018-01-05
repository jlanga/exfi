#!/usr/bin/env python3

from exfi.io import \
    index_fasta

import pandas as pd
import networkx as nx

bed6_cols = ["chrom", "start", "end", "name", "score", "strand"]

bed3records_empty = []
bed3records_simple = [
    ("ENSDART00000161035.1", 0, 326)
]
bed3records_complex = [
    ("ENSDART00000161035.1", 397, 472,),
    ("ENSDART00000165342.1", 1176, 1324),
    ("ENSDART00000161035.1", 0, 326),
    ("ENSDART00000165342.1", 125, 304),
    ("ENSDART00000165342.1", 746, 851),
    ("ENSDART00000165342.1", 974, 1097),
    ("ENSDART00000165342.1", 854, 886),
    ("ENSDART00000165342.1", 1098, 1175),
    ("ENSDART00000165342.1", 5, 127),
    ("ENSDART00000165342.1", 645, 746),
    ("ENSDART00000165342.1", 317, 460),
    ("ENSDART00000165342.1", 591, 650),
    ("ENSDART00000165342.1", 459, 592),
    ("ENSDART00000165342.1", 899, 953),
    ("ENSDART00000161035.1", 477, 523)
]

bed6df_empty = pd.DataFrame(
    columns=bed6_cols
)
bed6df_simple = pd.DataFrame(
    data=[("ENSDART00000161035.1", 0, 326, "ENSDART00000161035.1:0-326", 0, "+" )],
    columns=bed6_cols
)
bed6df_complex = pd.DataFrame(
        data=[
            ["ENSDART00000161035.1", 397, 472, "ENSDART00000161035.1:397-472", 0, "+"],
            ["ENSDART00000165342.1", 1176, 1324, "ENSDART00000165342.1:1176-1324", 0, "+"],
            ["ENSDART00000161035.1", 0, 326, "ENSDART00000161035.1:0-326", 0, "+"],
            ["ENSDART00000165342.1", 125, 304, "ENSDART00000165342.1:125-304", 0, "+"],
            ["ENSDART00000165342.1", 746, 851, "ENSDART00000165342.1:746-851", 0, "+"],
            ["ENSDART00000165342.1", 974, 1097, "ENSDART00000165342.1:974-1097", 0, "+"],
            ["ENSDART00000165342.1", 854, 886, "ENSDART00000165342.1:854-886", 0, "+"],
            ["ENSDART00000165342.1", 1098, 1175, "ENSDART00000165342.1:1098-1175", 0, "+"],
            ["ENSDART00000165342.1", 5, 127, "ENSDART00000165342.1:5-127", 0, "+"],
            ["ENSDART00000165342.1", 645, 746, "ENSDART00000165342.1:645-746", 0, "+"],
            ["ENSDART00000165342.1", 317, 460, "ENSDART00000165342.1:317-460", 0, "+"],
            ["ENSDART00000165342.1", 591, 650, "ENSDART00000165342.1:591-650", 0, "+"],
            ["ENSDART00000165342.1", 459, 592, "ENSDART00000165342.1:459-592", 0, "+"],
            ["ENSDART00000165342.1", 899, 953, "ENSDART00000165342.1:899-953", 0, "+"],
            ["ENSDART00000161035.1", 477, 523, "ENSDART00000161035.1:477-523", 0, "+"],
        ],
        columns=bed6_cols
    )\
    .sort_values(bed6_cols[0:3])

node2coords_empty = {}
node2coords_simple = {
    'ENSDART00000161035.1:0-326': ('ENSDART00000161035.1', 0, 326)
}
node2coords_complex = {
    'ENSDART00000161035.1:0-326': ('ENSDART00000161035.1', 0, 326),
    'ENSDART00000161035.1:397-472': ('ENSDART00000161035.1', 397, 472),
    'ENSDART00000161035.1:477-523': ('ENSDART00000161035.1', 477, 523),
    'ENSDART00000165342.1:5-127': ('ENSDART00000165342.1', 5, 127),
    'ENSDART00000165342.1:125-304': ('ENSDART00000165342.1', 125, 304),
    'ENSDART00000165342.1:317-460': ('ENSDART00000165342.1', 317, 460),
    'ENSDART00000165342.1:459-592': ('ENSDART00000165342.1', 459, 592),
    'ENSDART00000165342.1:591-650': ('ENSDART00000165342.1', 591, 650),
    'ENSDART00000165342.1:645-746': ('ENSDART00000165342.1', 645, 746),
    'ENSDART00000165342.1:746-851': ('ENSDART00000165342.1', 746, 851),
    'ENSDART00000165342.1:854-886': ('ENSDART00000165342.1', 854, 886),
    'ENSDART00000165342.1:899-953': ('ENSDART00000165342.1', 899, 953),
    'ENSDART00000165342.1:974-1097': ('ENSDART00000165342.1', 974, 1097),
    'ENSDART00000165342.1:1098-1175': ('ENSDART00000165342.1', 1098, 1175),
    'ENSDART00000165342.1:1176-1324': ('ENSDART00000165342.1', 1176, 1324)
}

path_empty = {}
path_simple = {"ENSDART00000161035.1": ("ENSDART00000161035.1:0-326",)}
path_complex = {
    "ENSDART00000161035.1": (
        "ENSDART00000161035.1:0-326",
        "ENSDART00000161035.1:397-472",
        "ENSDART00000161035.1:477-523",
    ),
    "ENSDART00000165342.1": (
        "ENSDART00000165342.1:5-127",
        "ENSDART00000165342.1:125-304",
        "ENSDART00000165342.1:317-460",
        "ENSDART00000165342.1:459-592",
        "ENSDART00000165342.1:591-650",
        "ENSDART00000165342.1:645-746",
        "ENSDART00000165342.1:746-851",
        "ENSDART00000165342.1:854-886",
        "ENSDART00000165342.1:899-953",
        "ENSDART00000165342.1:974-1097",
        "ENSDART00000165342.1:1098-1175",
        "ENSDART00000165342.1:1176-1324",
    )
}

overlaps_empty = {}
overlaps_simple = {}
overlaps_complex = {
    ("ENSDART00000161035.1:0-326", "ENSDART00000161035.1:397-472"): -71,
    ("ENSDART00000161035.1:397-472", "ENSDART00000161035.1:477-523"): -5,
    ("ENSDART00000165342.1:5-127", "ENSDART00000165342.1:125-304"): 2,
    ("ENSDART00000165342.1:125-304", "ENSDART00000165342.1:317-460"): -13,
    ("ENSDART00000165342.1:317-460", "ENSDART00000165342.1:459-592"): 1,
    ("ENSDART00000165342.1:459-592", "ENSDART00000165342.1:591-650"): 1,
    ("ENSDART00000165342.1:591-650", "ENSDART00000165342.1:645-746"): 5,
    ("ENSDART00000165342.1:645-746", "ENSDART00000165342.1:746-851"): 0,
    ("ENSDART00000165342.1:746-851", "ENSDART00000165342.1:854-886"): -3,
    ("ENSDART00000165342.1:854-886", "ENSDART00000165342.1:899-953"): -13,
    ("ENSDART00000165342.1:899-953", "ENSDART00000165342.1:974-1097"): -21,
    ("ENSDART00000165342.1:974-1097", "ENSDART00000165342.1:1098-1175"): -1,
    ("ENSDART00000165342.1:1098-1175", "ENSDART00000165342.1:1176-1324"): -1
}

splice_graph_empty = nx.DiGraph()

splice_graph_simple = nx.DiGraph()
splice_graph_simple.add_nodes_from(bed6df_simple["name"].tolist())
nx.set_node_attributes(
    G=splice_graph_simple,
    name="coordinates",
    values=node2coords_simple
)
for path in path_simple.values():
    splice_graph_simple.add_path(path)
nx.set_edge_attributes(
    G=splice_graph_simple,
    name="overlaps",
    values=overlaps_simple
)

splice_graph_complex = nx.DiGraph()
splice_graph_complex.add_nodes_from(bed6df_complex["name"].tolist())
nx.set_node_attributes(
    G=splice_graph_complex,
    name="coordinates",
    values=node2coords_complex
)
for path in path_complex.values():
    splice_graph_complex.add_path(path)
nx.set_edge_attributes(
    G=splice_graph_complex, name="overlaps", values=overlaps_complex
)


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
