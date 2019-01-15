#!/usr/bin/env python3

"""exfi.io.gfa1_to_exons.py: submodule to read a gfa1, extract the exons and
store it in fasta format"""

import pandas as pd

def gfa1_to_exons(fasta_out, gfa1_in):
    """Extract the exons in Fasta format"""
    with open(gfa1_in, "r") as gfa, open(fasta_out, "w") as fasta:
        segments = pd.DataFrame(
            data=[
                x.strip().split("\t")[0:3]
                for x in gfa.readlines() if x[0] == "S"
            ],
            columns=["RecordType", "Name", "Sequence"],
        )

        if segments.shape[0] == 0:
            return

        segments["fasta"] = ">" + segments["Name"] + "\n" + segments["Sequence"]
        segments.fasta.values.tofile(fasta, sep="\n", format="%s")
        fasta.write("\n")  # Final end line


def gfa1_to_gapped_transcripts(fasta_out, gfa1_in, gap_size=100):
    """Convert a GFA1 file to a gapped transcript file"""

    with open(gfa1_in, "r") as gfa, open(fasta_out, "w") as fasta:

        separator = gap_size * 'N'

        # Read only segments and paths
        data = [
            x.strip().split("\t")
            for x in gfa.readlines() if x[0] in set(["S", "P"])
        ]

        if not data:
            return

        # Create {node_id: nucleotide}
        node2sequence = pd.DataFrame(
            data=[x[0:3] for x in data if x[0] == "S"],
            columns=["RecordType", "Name", "Sequence"],
            )\
            .drop(columns="RecordType")\
            .set_index("Name")\
            .to_dict()["Sequence"]

        # Get the path info
        paths = pd.DataFrame(
            data=[x[0:4] for x in data if x[0] == "P"],
            columns=["RecordType", "PathName", "SegmentNames", "Overlaps"]
            )\
            .drop(columns=["RecordType", "Overlaps"])
        del data

        paths["SegmentNames"] = paths["SegmentNames"].str.replace("+", "")

        # Compose the sequence
        paths["gapped_sequence"] = paths\
            .SegmentNames\
            .str.split(',')\
            .map(lambda x: separator.join([node2sequence[y] for y in x]))
        del node2sequence

        # Create the fasta line
        paths["fasta"] = \
            ">" + paths.PathName + " " + paths.SegmentNames + "\n" + \
            paths.gapped_sequence

        # Dump everything
        paths.fasta.values.tofile(fasta, sep="\n", format="%s")
        fasta.write("\n")  # Final end line
