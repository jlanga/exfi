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
    # Get overlap and sequence data
    edge2overlap = nx.get_edge_attributes(G=splice_graph, name="overlap")
    node2seq = nx.get_node_attributes(G=splice_graph, name='sequence')

    # Prepare fasta for sealer
    # Make temporary fasta where to write sequences for sealer
    sealer_input = mkstemp()
    sequences_to_seal = list()
    for (node1, node2), overlap in edge2overlap.items():
        overlap = edge2overlap[(node1, node2)]
        if overlap < 0 and overlap < args["max_gap_size"]:
            identifier = node1 + "~" + node2
            sequence = Seq.Seq(node2seq[node1] + "N" * 100 + node2seq[node2])
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
    for corrected in SeqIO.parse(
        format="fasta", handle=handle
    ):
        node1, node2 = corrected.id.rsplit("_")[0].rsplit("_")[0].split("~")
        edge2fill[node1] = node2

    return edge2fill



def _sculpt_graph(splice_graph, edge2fill_tmp, transcriptome_index):
    """(nx.DiGraph, dict) -> nx.DiGraph

    Copy splice_graph, merge the nodes as dictated in edge2fill, return the
    sealed graph.
    """

    edge2fill = edge2fill_tmp

    while len(edge2fill) > 0:

        # Get involved nodes
        u = list(edge2fill.keys())[0]
        v = edge2fill[u]

        # Make new node name
        transcript1, start1, _ = _coordinates_to_variables(u)
        _, _, end2 = _coordinates_to_variables(v)
        new_node = transcript1 + ":" + str(start1) + "-" + str(end2)
        print(new_node)

        # Merge nodes
        splice_graph = nx.contracted_nodes( # It is collapsed into u!!
            G=splice_graph, u=u, v=v,
            self_loops=False
        )

        # Rename u
        splice_graph = nx.relabel_nodes(
            G=splice_graph,
            mapping={u: new_node},
            copy=False
        )

        # Update node2seq
        splice_graph.node[new_node]["sequence"] = str(transcriptome_index[transcript1][start1:end2].seq)
        # Update node2coord
        splice_graph.node[new_node]["coordinates"] = (transcript1, start1, end2)

        # Update edge2overlap
        #print(nx.get_edge_attributes(splice_graph, "overlap"))

        # Delete u from edge2fill
        del edge2fill[u]
        if v in edge2fill: # Update v if necessary
            edge2fill[new_node] = edge2fill[v]
            del edge2fill[v]



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
    return _sculpt_graph(splice_graph, edge2fill, transcriptome_index)
