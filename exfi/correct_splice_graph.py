#!/usr/bin/env python3

import networkx as nx

from subprocess import \
    Popen

from tempfile import \
    mkstemp

from os import remove

from Bio import \
    SeqIO, \
    Seq, \
    SeqRecord

from natsort import \
    natsorted




def _coordinates_to_variables(coordinates):
    """(string) -> (string, int, int)

    Convert a genomic coordinate of the form "TR_ID:start-end" into the tuple
    ("TR_ID", start, end).
    """
    transcript, start_end = coordinates.rsplit(":", 1)
    start_end = start_end.rsplit("-", 1)
    start = int(start_end[0])
    end = int(start_end[1])
    return transcript, start, end


def _get_node2sequence(splice_graph, transcriptome_dict):
    """(nx.DiGraph, dict) -> dict of str

    From the splice graph and a transcriptome, get the exon: sequence dictionary
    """
    node2coordinates = nx.get_node_attributes(
        G=splice_graph,
        name="coordinates"
    )

    node2sequence = {key: None for key in node2coordinates.keys()}

    for node, coordinate in node2coordinates.items():
        transcript_id, start, end = coordinate
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
        "bloom_filter": str,
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
        '--input-bloom', args["bloom_filter"],
        '--mask',
        '--output-prefix', sealer_output_prefix[1],
        '--verbose'
    ]

    # Execute
    p_sealer = Popen(c_sealer)
    p_sealer.communicate()

    # Clean files
    # remove(sealer_input_fn)
    # remove(sealer_output_prefix[1])
    remove(sealer_output_prefix[1] + "_log.txt")
    remove(sealer_output_prefix[1] + "_scaffold.fa")

    return sealer_output_prefix[1] + "_merged.fa"



def _collect_sealer_results(handle):
    """(str) -> dict

    Process extensions from sealer and return the computed extensions if a dict
    edge2fill = {(node1, node2): "acgt"}
    """
    # Collect results
    edge2fill = {}
    for corrected in SeqIO.parse(format="fasta", handle=handle):
        node1, node2 = corrected.id.rsplit("_")[0].rsplit("_")[0].split("~")
        edge2fill[node1] = node2

    return edge2fill



def _sculpt_graph(splice_graph, edge2fill):
    """(nx.DiGraph, dict) -> nx.DiGraph

    Copy splice_graph, merge the nodes as dictated in edge2fill, return the
    sealed graph.
    """

    while len(edge2fill) > 0:

        # Get nodes to modify
        u = natsorted(edge2fill.keys())[0]
        v = edge2fill[u]

        # Compose new names and coordinates
        u_transcript, u_start, u_end = splice_graph.node[u]["coordinates"]
        v_transcript, v_start, v_end = splice_graph.node[v]["coordinates"]
        n_coordinates = (u_transcript, u_start, v_end)
        n = "{0}:{1}-{2}".format(*n_coordinates)

        # Insert new node
        splice_graph.add_node(n)
        splice_graph.node[n]["coordinates"] = n_coordinates

        # link pred(u) to n, and overlaps
        for predecessor in splice_graph.predecessors(u):
            splice_graph.add_edge(u=predecessor, v=n)
            splice_graph[predecessor][n]['overlaps'] = splice_graph[predecessor][u]['overlaps']

        # Attach n to succ(v)
        for successor in splice_graph.successors(v):
            splice_graph.add_edge(u=n, v=successor)
            splice_graph[n][successor]['overlaps'] = splice_graph[v][successor]['overlaps']

        # Delete u, delete v
        splice_graph.remove_nodes_from(nodes=[u, v])

        # Update dict of edges to fill
        if v in edge2fill:  # Update v if necessary
            edge2fill[n] = edge2fill[v]
            del edge2fill[v]
        del edge2fill[u]

    return splice_graph



def correct_splice_graph(splice_graph, args):
    """(nx.DiGraph, str, int, int) -> nx.DiGraph

    Try to correct small gaps (SNPs and indels) with abyss-sealer

    args = {
        "kmer": int,
        "max_gap_size": int,
        "bloom_filter": str,
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

    transcriptome_index = SeqIO.index(
        filename=args["input_fasta"],
        format="fasta"
    )

    # Compute the sealed splice graph
    return _sculpt_graph(splice_graph, edge2fill)
