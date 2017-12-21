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
    transcript, start_end = coordinates.rsplit(":")
    start_end = start_end.split("-")
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
    remove(sealer_input_fn)
    remove(sealer_output_prefix[1])
    remove(sealer_output_prefix[1] + "_log.txt")
    remove(sealer_output_prefix[1] + "_scaffold.fa")

    return sealer_output_prefix[1] + "_merged.fa"



def _collect_sealer_results(handle):
    """(str) -> dict

    Process extensions from sealer and return the computed extensions if a dict
    edge2fill = {(node1, node2): "acgt"}
    """
    # Collect results
    # translate ambiguous nucleotides into the lexicografically smallest one
    iupac_translator = { # Alphabetically
        'a': 'a', 'c': 'c', 'g': 'g', 't': 't',
        'r': 'r', 'y': 'c', 's': 'c', 'w': 'a', 'k': 'g', 'm': 'a',
        'b': 'c', 'd': 'a', 'h': 'v',
        'n': 'a'
    }

    # Collect results
    edge2fill = {}
    for corrected in SeqIO.parse(
        format="fasta", handle=handle
    ):
        node1, node2 = corrected.id.rsplit("_")[0].rsplit("_")[0].split("~")
        edge2fill[(node1, node2)] =  ''.join(
            [iupac_translator[x] for x in str(corrected.seq) if x.islower()]
        )

    # Clean sealer inputs and outputs
    remove(handle)

    return edge2fill



def _sculpt_graph(splice_graph, edge2fill):
    """(nx.DiGraph, dict) -> nx.DiGraph

    Copy splice_graph, merge the nodes as dictated in edge2fill, return the
    sealed graph.
    """
    # Copy old graph
    sealed_splicegraph = splice_graph.copy()

    # Extract node data
    node2seq = nx.get_node_attributes(G=sealed_splicegraph, name="sequence")
    node2coord = nx.get_node_attributes(G=sealed_splicegraph, name="coordinates")

    # Extract edge data
    edge2overlap = nx.get_edge_attributes(G=sealed_splicegraph, name="overlap")

    # Modify sealed nodes
    for (node1, node2), fill in edge2fill.items():

        # Manipulate nodes:
        # Get old node data
        transcript, start, _ = _coordinates_to_variables(node1)
        _, _, end = _coordinates_to_variables(node2)

        # Add new node name: transcript1:start1-end2
        new_node = "{transcript}:{start}-{end}".format(
            transcript=transcript,
            start=start,
            end=end
        )
        sealed_splicegraph.add_node(new_node)

        # Compose the new sequence
        new_sequence = node2seq[node1] + fill.upper() + node2seq[node2]
        node2seq[new_node] = new_sequence

        # Add new coordinates
        node2coord[new_node] = tuple([transcript, start, end])

        # Manipulate edges:
        predecessor_list = list(sealed_splicegraph.predecessors(n=node1))
        successor_list = list(sealed_splicegraph.successors(n=node2))

        # Test node1 is first node in transcript: no predecessor
        if predecessor_list:
            predecessor = predecessor_list[0]
            # Add new edge
            sealed_splicegraph.add_edge(u=predecessor, v=new_node)
            # Get andc modify overlap data
            overlap = edge2overlap[(predecessor, node1)]
            edge2overlap[(predecessor, new_node)] = overlap

        # Test node2 is last node in transcript: No successor
        if successor_list:
            successor = successor_list[0]
            # Add new, delete old
            sealed_splicegraph.add_edge(u=new_node, v=successor)
            # Get overlap data
            overlap = edge2overlap[(node2, successor)]
            edge2overlap[(new_node, successor)] = overlap

    # clean up
    nodes_to_clean = set()
    edges_to_clean = set()

    for (node1, node2) in edge2fill.keys():
        nodes_to_clean.add(node1)
        nodes_to_clean.add(node2)
        edges_to_clean.add((node1, node2))

    nodes_to_clean = list(nodes_to_clean)
    edges_to_clean = list(edges_to_clean)

    # Clean node data
    sealed_splicegraph.remove_nodes_from(nodes_to_clean)
    sealed_splicegraph.remove_edges_from(edges_to_clean)
    for node in nodes_to_clean:
        del node2seq[node]
        del node2coord[node]
    for edge in edges_to_clean:
        del edge2overlap[edge]


    # overwrite node2coord, node2seq, edge2overlap
    nx.set_node_attributes(G=sealed_splicegraph, name="coordinates", values=node2coord)
    nx.set_node_attributes(G=sealed_splicegraph, name="sequence", values=node2seq)
    nx.set_edge_attributes(G=sealed_splicegraph, name="overlap", values=edge2overlap)

    return sealed_splicegraph


def correct_splice_graph(splice_graph, args):
    """(nx.DiGraph, str, int, int) -> nx.DiGraph

    Try to correct small gaps (SNPs and indels) with abyss-sealer

    args = {
        "kmer": int,
        "max_gap_size": int,
        "bloom_filter": str,
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

    # Compute the sealed splice graph
    return _sculpt_graph(splice_graph, edge2fill)
