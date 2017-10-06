#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import networkx as nx
import pandas as pd

from exfi.build_splicegraph import \
    exons_to_df, \
    transcript_to_path


def write_gfa1(splice_graph, transcript_index, exons, filename):
    """(dict_of_seqrecords, dict_of_seqrecords, str) -> None

    Write splice graph to filename in GFA 1 format
    """
    with open(filename, "w") as gfa:

        # Header
        gfa.write("H\tVN:Z:1.0\n")

        # Nodes
        ## Transcripts (commented)
        for transcript_id in sorted(transcript_index.keys()):
            transcript = transcript_index[transcript_id]
            gfa.write(
                "#S\t{identifier}\t{sequence}\tLN:i:{length}\n".format(
                    identifier=transcript.id,
                    sequence=str(transcript.seq),
                    length=len(transcript.seq)
                )
            )

        ## Exons
        node2seq = nx.get_node_attributes(
            G=splice_graph,
            name='sequence'
        )
        for node in sorted(splice_graph.nodes()):
            gfa.write(
                "S\t{node}\t{sequence}\tLN:i:{length}\n".format(
                    node=node,
                    sequence=node2seq[node],
                    length=len(node2seq[node])
            )
        )

        # Add coordinates as containments
        exon2coordinates = nx.get_node_attributes(
            G=splice_graph,
            name='coordinates'
        )
        for exon_id in sorted(exon2coordinates.keys()):
            coordinates = exon2coordinates[exon_id]
            for coordinate in coordinates:
                transcript_id, start, end = coordinate
                gfa.write(
                    "C\t{container_id}\t{container_orient}\t"
                    "{contained_id}\t{contained_orient}\t"
                    "{position}\t{overlap}\n".format(
                        container_id=transcript_id,
                        container_orient="+",
                        contained_id=exon_id,
                        contained_orient="+",
                        position=start,
                        overlap=str(len(node2seq[exon_id])) + "M"
                    )
                )

        # Edges
        edge2overlap = nx.get_edge_attributes(
            G=splice_graph,
            name="overlap"
        )
        for edge in sorted(splice_graph.edges()):
            node1, node2 = edge
            overlap = edge2overlap[edge]
            # is an overlap or a gap
            if overlap >= 0:
                overlap = "{}M".format(overlap)
            else:
                overlap = "{}G".format(-overlap)
            gfa.write(
                "L\t{node1}\t{orientation1}\t{node2}\t{orientation2}\t{overlap}\n".format(
                    node1=node1,
                    orientation1="+",
                    node2=node2,
                    orientation2="+",
                    overlap = overlap
                )
            )

        # Paths
        exon_df = exons_to_df(exons)
        paths = transcript_to_path(exon_df)
        for transcript_id in sorted(paths.keys()):
            gfa.write(
                "P\t{transcript_id}\t{path}\n".format(
                    transcript_id=transcript_id,
                    path=",".join([exon+"+" for exon in paths[transcript_id]])
                )
            )


def gfa1_to_exons(gfa_in_fn, fasta_out_fn, soft_mask_overlaps=False):

    with open(gfa_in_fn, "r") as gfa_in:

        # Get the seqrecord
        exon_dict = {}
        coordinate_dict = {}
        overlap_dict = {}

        for line in gfa_in:
            line = line.strip()

            # Nodes
            if line.startswith("S"):
                _, exon_id, sequence, _ = line.split("\t")
                exon_dict[exon_id] = SeqRecord(id=exon_id, description="", seq=Seq(sequence))

            # Containments (coordinates)
            if line.startswith("C"):
                _, container_id, container_orientation, contained_id, contained_orientation, position, overlap = line.split("\t")
                if contained_id not in coordinate_dict:
                    coordinate_dict[contained_id] = []
                coordinate_dict[contained_id].append(
                    "{transcript_id}:{start}-{end}".format(
                        transcript_id = container_id,
                        start=position,
                        end=int(position) + int(overlap[:-1])
                    )
                )

            # Overlaps
            if line.startswith("L") and soft_mask_overlaps:
                _, start, _, end, _, overlap = line.split("\t")
                letter = overlap[-1]
                overlap = int(overlap[:-1])
                if overlap > 0 and letter == "M":
                    # node1
                    start_seq = exon_dict[start]
                    start_seq = start_seq[:-overlap] + start_seq[-overlap:].lower()
                    exon_dict[start] = start_seq

                    # node2
                    end_seq = exon_dict[end]
                    end_seq = end_seq[:overlap].lower() + end_seq[overlap:]
                    exon_dict[end] = end_seq

        # Join sequence and coordinates
        for exon_id, exon_record in exon_dict.items():
            exon_record.description = " ".join(coordinate_dict[exon_id])

        # Write to fasta
        SeqIO.write(format="fasta", handle=fasta_out_fn, sequences=exon_dict.values())



def gfa1_to_gapped_transcript(gfa_in, fasta_out, number_of_ns=100, soft_mask_overlaps=False):

    # Process gfa
    with open(gfa_in, "r") as gfa_in:
        exon_dict = {} # exon_id: sequence as str
        path_dict = {} # transcript_id : [exon_id1, exon_idN]
        overlap_dict = {} # (start, end): int

        for line in gfa_in:

            line = line.strip()

            # Process Segments
            if line.startswith("S"):
                _, exon_id, sequence = line.split("\t")[0:3]
                exon_dict[exon_id] = sequence

            # Process links
            if line.startswith("L") and soft_mask_overlaps:
                _, start, _, end, _, overlap = line.split("\t")
                letter = overlap[-1]
                overlap = int(overlap[:-1])
                if letter == "M" and overlap > 0:
                    overlap_dict[(start, end)] = overlap

            # Process paths
            if line.startswith("P"):
                _, path_id, path = line.split("\t")[0:3]
                path_list = path.split(",")
                path_list = [exon_id[0:-1] for exon_id in path_list]
                path_dict[path_id] = path_list

    # Mask all overlaps in the exon_dict
    if soft_mask_overlaps:
        for edge, overlap in overlap_dict.items():
            start, end = edge
            # node1
            start_seq = exon_dict[start]
            start_seq = start_seq[:-overlap] + start_seq[-overlap:].lower()
            exon_dict[start] = start_seq

            # node2
            end_seq = exon_dict[end]
            end_seq = end_seq[:overlap].lower() + end_seq[overlap:]
            exon_dict[end] = end_seq

    for k, v in exon_dict.items():
        print(k, ":", v)

    # Compose path as seq
    ns = "N" * number_of_ns
    records = []
    for transcript_id, exon_list in path_dict.items():
        exon_seqs = [exon_dict[exon_id] for exon_id in exon_list]
        records.append(SeqRecord(
            id=transcript_id,
            description=",".join(exon_list),
            seq=Seq(ns.join(exon_seqs))
        ))

    # Write
    SeqIO.write(
        format="fasta",
        sequences=records,
        handle=fasta_out
    )
