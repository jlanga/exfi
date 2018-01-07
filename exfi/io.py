#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import networkx as nx

from exfi.build_splice_graph import \
    bed3_records_to_bed6df, \
    bed6df_to_path2node

from itertools import chain

def _clean_seqrecord(seqrecord):
    """Delete the identifier from the description"""
    seqrecord.description = " ".join(seqrecord.description.split(" ")[1:])
    return seqrecord



def _clean_index(index):
    """Clean all elements from an indexed fasta"""
    index_clean = {}
    for key, value in index.items():
        index_clean[key] = _clean_seqrecord(value)
    return index_clean



def index_fasta(filename):
    """Create a fasta dict, with clean descriptions, key=id, value=seqrecord"""
    index = SeqIO.index(
        filename=filename,
        format="fasta"
    )
    return _clean_index(index)



def _segments_to_exon_dict(segments):
    """(list of lists) -> dict

    Conver a list of ["S", "ident", "seq", *whatever] to a dict {ident: SeqRecord}
    """
    exon_dict = {}
    for segment in segments:
        exon_id, sequence = segment[1:3]
        exon_dict[exon_id] = SeqRecord(
            id=exon_id,
            description="",
            seq=Seq(sequence)
        )
    return exon_dict



def _containments_to_coordinate_dict(containments):
    """(list of lists) -> dict

    Convert a list of ["C", transcript_id, +, exon_id, +, position, overlap ] to a dict of
    {exon_id: [coordinates wrt transcriptome]}
    """
    coordinate_dict = {}  # {exon: [(tr_id1:start-end),..., (id2:start-end)]}

    for containment in containments:

        # Take what is needed
        _, transcript_id, _, exon_id, _, position, overlap = containment
        if exon_id not in coordinate_dict.keys():
            coordinate_dict[exon_id] = []

        coordinate_dict[exon_id].append(
            "{transcript_id}:{start}-{end}".format(
                transcript_id=transcript_id,
                start=position,
                end=int(position) + int(overlap[:-1])
            )
        )

    return coordinate_dict



def _links_to_overlap_dict(links):
    """(list of lists) -> dict

    Convert a list of ["L", transcript_id, +, exon_id, +, position, overlap ] to a dict of
    {exon_id: [coordinates wrt transcriptome]}
    """
    overlap_dict = {}

    for link in links:
        _, start, _, end, _, overlap = link
        overlap_dict[(start, end)] = overlap
    return overlap_dict



def _paths_to_path_dict(paths):
    """(list_of_lists) -> dict

    Convert a list of ["P", transcript_id, [exon], overlap ] to a dict of
    {exon_id: [coordinates wrt transcriptome]}
    """
    path_dict = {}

    for path in paths:
        _, path_id, path, *_ = path
        path_list = path.split(",")
        path_list = [exon_id[0:-1] for exon_id in path_list]  # Delete the +s
        path_dict[path_id] = path_list

    return path_dict



def read_gfa1(filename):

    with open(filename, "r") as gfain:

        header = []
        segments = []
        links = []
        containments = []
        paths = []


        for line in gfain:

            line = line.strip().split("\t")

            if line[0] == "H": header.append(line)
            elif line[0] == "S": segments.append(line)
            elif line[0] == "L": links.append(line)
            elif line[0] == "C": containments.append(line)
            elif line[0] == "P": paths.append(line)

    return {
        "header": header,
        "segments": segments,
        "links": links,
        "containments": containments,
        "paths": paths
    }



def _process_overlap_cigar(cigar_string):
    """Process a simple CIGAR string (number, letter)"""
    return [cigar_string[-1], int(cigar_string[:-1])]



def _soft_mask_right(string, n):
    """Soft mask the rightmost n bases"""
    return string[:-n] + string[-n:].lower()



def _soft_mask_left(string, n):
    """Soft mask the leftmost n bases"""
    return string[:n].lower() + string[n:]



def _soft_mask(exon_dict, overlap_dict):
    """Soft mask all overlaps in the exon_dict"""
    for (start, end), overlap in overlap_dict.items():
        letter, overlap = _process_overlap_cigar(overlap)
        if letter == "M" and overlap > 0:
            exon_dict[start] = _soft_mask_right(exon_dict[start], overlap)
            exon_dict[end] = _soft_mask_left(exon_dict[end], overlap)
    return exon_dict



def _hard_mask_right(string, n):
    """Hard mask the rightmost n bases"""
    return string[:-n] + "N" * n



def _hard_mask_left(string, n):
    """Hard mask the leftmost n bases"""
    return "N" * n + string[n:]



def _hard_mask(exon_dict, overlap_dict):
    """Hard mask all overlaps in the exon_dict"""
    for (start, end), overlap in overlap_dict.items():
        letter, overlap = _process_overlap_cigar(overlap)
        if letter == "M" and overlap > 0:
            exon_dict[start] = _hard_mask_right(exon_dict[start], overlap)
            exon_dict[end] = _hard_mask_left(exon_dict[end], overlap)
    return exon_dict



def _mask(exon_dict, overlap_dict, soft_mask_overlaps=False, hard_mask_overlaps=False):
    """If any of the soft mask or hard mask are activated, mask"""
    if soft_mask_overlaps == True and hard_mask_overlaps == True:
        raise Exception("I can't soft mask and hard mask at the same time, dude!")
    if soft_mask_overlaps: exon_dict = _soft_mask(exon_dict, overlap_dict)
    if hard_mask_overlaps: exon_dict = _hard_mask(exon_dict, overlap_dict)



def _compute_segments(splice_graph, transcriptome_dict):
    """(nx.DiGraph, dict) -> list

    Compute the segment lines:
    S node_id sequence length
    """
    node2coords = nx.get_node_attributes(G=splice_graph, name="coordinates")
    for node_id, coordinate in node2coords.items():
        node, start, end = coordinate
        sequence = str(transcriptome_dict[node].seq[start:end])
        yield "S\t{node}\t{sequence}\tLN:i:{length}\n".format(
            node=node_id,
            sequence=sequence,
            length=len(sequence)
        )



def _compute_links(splice_graph):
    """(nx.DiGraph) -> list

    Compute the link lines:
    L start orientation end orientation overlap
    """
    # Edges
    edge2overlap = nx.get_edge_attributes(
        G=splice_graph,
        name="overlaps"
    )

    for (node1, node2), overlap in edge2overlap.items():
        # is an overlap or a gap
        if overlap >= 0:
            overlap = "{}M".format(overlap)
        else:
            overlap = "{}G".format(-overlap)

        yield "L\t{node1}\t{orientation1}\t{node2}\t{orientation2}\t{overlap}\n".format(
            node1=node1,
            orientation1="+",
            node2=node2,
            orientation2="+",
            overlap = overlap
        )



def _compute_containments(splice_graph, transcriptome_dict):
    """(nx.DiGraph, dict of SeqRecords) -> list

    Compute the containment lines (w.r.t. transcriptome):
    C container orientation contained orientation position overlap
    """
    # Extract from the graph necessary data
    node2coordinates = nx.get_node_attributes(
        G=splice_graph,
        name='coordinates'
    )
    for node, coordinates in node2coordinates.items():
        transcript_id, start, end = coordinates
        sequence = str(transcriptome_dict[transcript_id].seq[start:end])
        yield "C\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
            transcript_id, "+",
            node, "+",
            start, str(len(sequence)) + "M"
        )



def _compute_paths(splice_graph):
    """(nx.DiGraph) -> lists

    Compute the paths in the splice graph:
    P transcript_id [node1, ..., nodeN]
    """
    bed3records = nx.get_node_attributes(
        G=splice_graph,
        name="coordinates"
    ).values()

    bed6df = bed3_records_to_bed6df(bed3records)
    path2nodes = bed6df_to_path2node(bed6df)
    for transcript_id, path in path2nodes.items():
        yield "P\t{transcript_id}\t{path}\n".format(
            transcript_id=transcript_id,
            path=",".join([node + "+" for node in path])
        )




def write_gfa1(splice_graph, transcriptome_dict, filename):
    """(dict_of_seqrecords, dict_of_seqrecords, str) -> None

    Write splice graph to filename in GFA 1 format
    """
    header = ["H\tVN:Z:1.0\n"]

    segments = list(_compute_segments(splice_graph, transcriptome_dict))
    links = list(_compute_links(splice_graph))
    containments = list(_compute_containments(splice_graph, transcriptome_dict))
    paths = list(_compute_paths(splice_graph))
    with open(filename, "w") as gfa1_out:
        gfa1_out.writelines(chain(
            header, segments, links, containments, paths
        ))


def gfa1_to_exons(gfa_in_fn, fasta_out_fn, soft_mask_overlaps=False, hard_mask_overlaps=False):

    gfa1 = read_gfa1(gfa_in_fn)

    exon_dict = _segments_to_exon_dict(gfa1["segments"])
    overlap_dict = _links_to_overlap_dict(gfa1["links"])

    # Mask if necessary
    _mask(exon_dict, overlap_dict, soft_mask_overlaps, hard_mask_overlaps)

    # Add coordinate information to description
    for exon_id, exon_record in exon_dict.items():
        exon_record.description = " ".join(coordinate_dict[exon_id])

    # Write to fasta
    SeqIO.write(format="fasta", handle=fasta_out_fn, sequences=exon_dict.values())



def _compose_paths(exon_dict, path_dict, number_of_ns):
    ns = "N" * number_of_ns
    for transcript_id, exon_list in sorted(path_dict.items()):
        exon_seqs = [str(exon_dict[exon_id].seq) for exon_id in exon_list]
        yield SeqRecord(
            id=transcript_id,
            description=",".join(exon_list),
            seq=Seq(ns.join(exon_seqs))
        )


def gfa1_to_gapped_transcript(
    gfa_in, fasta_out, number_of_ns=100, soft_mask_overlaps=False, hard_mask_overlaps=False):

    if soft_mask_overlaps == True and hard_mask_overlaps == True:
        raise Exception("I can't soft mask and hard mask at the same time, dude!")

    # Process
    gfa1 = read_gfa1(gfa_in)
    exon_dict = _segments_to_exon_dict(gfa1["segments"])
    coordinate_dict = _containments_to_coordinate_dict(gfa1["containments"])
    overlap_dict = _links_to_overlap_dict(gfa1["links"])
    path_dict = _paths_to_path_dict(gfa1["paths"])


    # Mask if necessary
    _mask(exon_dict, overlap_dict, soft_mask_overlaps, hard_mask_overlaps)

    composed_paths = _compose_paths(exon_dict, path_dict, number_of_ns)

    # Write
    SeqIO.write(format="fasta", sequences=composed_paths, handle=fasta_out)
