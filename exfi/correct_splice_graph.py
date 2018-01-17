#!/usr/bin/env python3

"""
exfi.correct_splice_graph.py: functions to take a splice graph and try to fill hypothetical gaps
with abyss-sealer.
"""

from subprocess import \
    Popen

from tempfile import \
    mkstemp

from os import remove

import networkx as nx

from Bio import \
    SeqIO, \
    Seq, \
    SeqRecord

from natsort import \
    natsorted



def _get_node2sequence(splice_graph, transcriptome_dict):
    """(nx.DiGraph, dict) -> dict of str

    From the splice graph and a transcriptome, get the exon: sequence dictionary
    """
    node2coordinates = nx.get_node_attributes(
        G=splice_graph,
        name="coordinates"
    )

    node2sequence = {key: None for key in node2coordinates.keys()}

    for node, coordinates in node2coordinates.items():
        for (transcript_id, start, end) in coordinates:
            sequence = str(transcriptome_dict[transcript_id][start:end].seq)
            node2sequence[node] = sequence
    return node2sequence


def _prepare_sealer(splice_graph, args):
    """ (nx.DiGraph, dict_of_parameters) -> str

    Prepare fasta file with candidates to be filled with sealer. Return the path
    of the fasta file to be sealed.

    args = {
        "kmer": int,
        "max_gap_size": int, <- only one used
        "input_bloom": str,
    }
    """
    transcriptome_dict = SeqIO.index(args["input_fasta"], format="fasta")
    # Get overlap and sequence data
    edge2overlap = nx.get_edge_attributes(G=splice_graph, name="overlaps")
    node2sequence = _get_node2sequence(
        splice_graph=splice_graph,
        transcriptome_dict=transcriptome_dict
    )

    # Prepare fasta for sealer
    # Make temporary fasta where to write sequences for sealer
    sealer_input = mkstemp()
    sequences_to_seal = list()
    for (node1, node2), overlap in edge2overlap.items():
        overlap = edge2overlap[(node1, node2)]
        if overlap < 0 and overlap < args["max_gap_size"]:
            identifier = node1 + "~" + node2
            sequence = Seq.Seq(node2sequence[node1] + "N" * 100 + node2sequence[node2])
            seqrecord = SeqRecord.SeqRecord(
                id=identifier,
                description="",
                seq=sequence
            )
            sequences_to_seal.append(seqrecord)

    SeqIO.write(
        format="fasta",
        handle=sealer_input[1],
        sequences=sequences_to_seal
    )

    return sealer_input[1]



def _run_sealer(sealer_input_fn, args):
    """(str, dict) -> str

    Run abyss-sealer with the parameters in args, and the scaffold in
    sealer_input.

    args = {
        "kmer": int,
        "max_gap_size": int,
        "input_bloom": str,
    }
    """
    # Run sealer
    sealer_output_prefix = mkstemp()
    c_sealer = [
        'abyss-sealer',
        '--input-scaffold', sealer_input_fn,
        '--flank-length', str(args["kmer"]),
        '--max-gap-length', str(args["max_gap_size"]),
        '--kmer', str(args["kmer"]),
        '--fix-errors',
        '--input-bloom', args["input_bloom"],
        '--mask',
        '--output-prefix', sealer_output_prefix[1],
        '--verbose'
    ]

    # Execute
    p_sealer = Popen(c_sealer)
    p_sealer.communicate()

    # Clean files
    remove(sealer_output_prefix[1] + "_log.txt")
    remove(sealer_output_prefix[1] + "_scaffold.fa")
    remove(sealer_output_prefix[1])
    return sealer_output_prefix[1] + "_merged.fa"



def _collect_sealer_results(handle):
    """(str) -> dict

    Process extensions from sealer and return the computed extensions if a dict
    edge2fill = {(node1, node2): "acgt"}
    """
    # Collect results
    edge2fill = {}
    for corrected in SeqIO.parse(format="fasta", handle=handle):
        node1, node2 = corrected.id.rsplit("_", 2)[0].split("~")
        edge2fill[node1] = node2

    return edge2fill



def _sculpt_graph(splice_graph, edge2fill):
    """(nx.DiGraph, dict) -> nx.DiGraph

    Copy splice_graph, merge the nodes as dictated in edge2fill, return the
    sealed graph.
    """

    while edge2fill:

        # Get nodes to modify
        node_u = natsorted(edge2fill.keys())[0]
        node_v = edge2fill[node_u]

        # Compose new names and coordinates
        u_transcript, u_start, _ = splice_graph.node[node_u]["coordinates"][0]
        _, _, v_end = splice_graph.node[node_v]["coordinates"][0]
        n_coordinates = (u_transcript, u_start, v_end)
        node_n = "{0}:{1}-{2}".format(*n_coordinates)

        # Insert new node
        splice_graph.add_node(node_n)
        splice_graph.node[node_n]["coordinates"] = (n_coordinates,)

        # link pred(u) to n, and overlaps
        for predecessor in splice_graph.predecessors(node_u):
            splice_graph.add_edge(u=predecessor, v=node_n)
            splice_graph[predecessor][node_n]['overlaps'] = \
                splice_graph[predecessor][node_u]['overlaps']

        # Attach n to succ(v)
        for successor in splice_graph.successors(node_v):
            splice_graph.add_edge(u=node_n, v=successor)
            splice_graph[node_n][successor]['overlaps'] = \
                splice_graph[node_v][successor]['overlaps']

        # Delete u, delete v
        splice_graph.remove_nodes_from(nodes=[node_u, node_v])

        # Update dict of edges to fill
        if node_v in edge2fill:  # Update v if necessary
            edge2fill[node_n] = edge2fill[node_v]
            del edge2fill[node_v]
        del edge2fill[node_u]

    return splice_graph



def correct_splice_graph(splice_graph, args):
    """(nx.DiGraph, str, int, int) -> nx.DiGraph

    Try to correct small gaps (SNPs and indels) with abyss-sealer

    args = {
        "kmer": int,
        "max_gap_size": int,
        "input_bloom": str,
        "input_fasta": str
    }
    """
    # Compose fasta with candidates to be filled
    sealer_input_fn = _prepare_sealer(
        splice_graph=splice_graph,
        args=args
    )

    # Run sealer
    sealer_output_fn = _run_sealer(
        sealer_input_fn=sealer_input_fn,
        args=args
    )

    # Collect sealer results
    edge2fill = _collect_sealer_results(handle=sealer_output_fn)

    remove(sealer_input_fn)
    remove(sealer_output_fn)

    # Compute the sealed splice graph
    return _sculpt_graph(splice_graph, edge2fill)
