#!/usr/bin/env python3

"""
exfi.correct_splice_graph.py: functions to take a splice graph and try to fill hypothetical gaps
with abyss-sealer.
"""

import logging

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

from exfi.io.fasta_to_dict import \
    fasta_to_dict

def _get_node2sequence(splice_graph: dict, transcriptome_dict: dict) -> dict:
    """From the splice graph and a transcriptome, get the exon: sequence dictionary"""
    logging.info("\tComputing the exon to sequence dictionary")
    node2sequence = {}

    node2coordinates = nx.get_node_attributes(
        G=splice_graph,
        name="coordinates"
    )

    for node, coordinates in node2coordinates.items():
        for (transcript_id, start, end) in coordinates:
            sequence = transcriptome_dict[transcript_id][start:end]
            node2sequence[node] = sequence
    return node2sequence



def _prepare_sealer(splice_graph_dict: nx.DiGraph, args: dict) -> str:
    """Prepare fasta file with candidates to be filled with sealer. Return the path
    of the fasta file to be sealed.

    Candidates to be sealed:
        - pairs of exons with mall gaps (size <= max_gap_size)
        - pairs of exons with any kind of positive overlap

    args = {
        "input_fasta": str,
        "kmer": int,
        "max_gap_size": int,
        "input_bloom": str,
        "max_fp_bases": int
    }
    """

    logging.info("\tPreparing input for abyss-sealer")
    transcriptome_dict = fasta_to_dict(args["input_fasta"])

    # Prepare fasta for sealer
    # Make temporary fasta where to write sequences for sealer
    sealer_input = mkstemp()
    sequences_to_seal = list()

    for splice_graph in splice_graph_dict.values():

        # Get overlap and sequence data
        edge2overlap = nx.get_edge_attributes(G=splice_graph, name="overlaps")
        node2sequence = _get_node2sequence(
            splice_graph=splice_graph,
            transcriptome_dict=transcriptome_dict
        )

        for (node1, node2), overlap in edge2overlap.items():
            overlap = edge2overlap[(node1, node2)]
            if overlap < 0 and overlap <= args["max_gap_size"]:  # Small gap
                identifier = node1 + "~" + node2
                sequence = Seq.Seq(
                    node2sequence[node1][0:-args["max_fp_bases"]] \
                    + "N" * 100 \
                    + node2sequence[node2][args["max_fp_bases"]:]
                )
                seqrecord = SeqRecord.SeqRecord(
                    id=identifier,
                    description="",
                    seq=sequence
                )
                sequences_to_seal.append(seqrecord)
            elif overlap >= 0:
                # Trim overlap bases from one of the threads
                identifier = node1 + "~" + node2
                # Cut overlap from one end, again from the other, and from both to see what happens

                ## Both
                sequence = Seq.Seq(
                    node2sequence[node1][0:-overlap-args["max_fp_bases"]] \
                    + "N" * 100 \
                    + node2sequence[node2][overlap+args["max_fp_bases"]:]
                )
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



def _run_sealer(sealer_input_fn: str, args: dict) -> str:
    """Run abyss-sealer with the parameters in args, and the scaffold in
    sealer_input.

    args = {
        "kmer": int,
        "max_gap_size": int,
        "input_bloom": str,
    }
    """
    logging.info("\tRunning abyss-sealer")
    # Run sealer
    sealer_output_prefix = mkstemp()
    c_sealer = [
        'abyss-sealer',
        '--input-scaffold', sealer_input_fn,
        '--flank-length', str(args["kmer"]),
        '--max-gap-length', "30",
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



def _collect_sealer_results(handle: str) -> set:
    """Process extensions from sealer and return the computed extensions in a set of pairs of
    tuples: set((node1, node2), ... , (nodeN-1, nodeN))"""
    logging.info("\tCollecting abyss-sealer results")
    # Collect results
    filled_edges = set()
    for corrected in fasta_to_dict(handle).keys():
        node1, node2 = corrected.rsplit("_", 2)[0].split("~")
        filled_edges.add((node1, node2))
    return filled_edges



# def _filled_edges_by_transcript(splice_graph: nx.DiGraph, filled_edges: str) -> dict:
#     """Split the edge2fill by the transcript they belong. Result is
#     dict(transcript_id: set)"""
#     logging.info("\tSplitting sealer results by transcript")
#     filled_edges_by_transcript = {}
#     for node_u, node_v in filled_edges:
#         transcript = splice_graph.nodes[node_u]["coordinates"][0][0]
#         if transcript not in filled_edges_by_transcript:
#             filled_edges_by_transcript[transcript] = set()
#         filled_edges_by_transcript[transcript].add((node_u, node_v))
#     return filled_edges_by_transcript



def _filled_edges_by_transcript(filled_edges: str) -> dict:
    """Split the edge2fill by the transcript they belong. Result is
    dict(transcript_id: set)"""
    logging.info("\tSplitting sealer results by transcript")
    filled_edges_by_transcript = {}
    for node_u, node_v in filled_edges:
        transcript = node_u.rsplit(":")[0]
        if transcript not in filled_edges_by_transcript:
            filled_edges_by_transcript[transcript] = set()
        filled_edges_by_transcript[transcript].add((node_u, node_v))
    return filled_edges_by_transcript



def _rename_nodes_from_collapse(quotient_graph: nx.DiGraph) -> dict:
    """Compose the new_node ids from nx.quotient to str or tuples of strs"""
    logging.info("\tRenaming collapsed nodes")
    # Main dict
    mapping = {  # Old -> New
        node_id: tuple(natsorted(node for node in node_id))
        for node_id in quotient_graph.nodes()
    }
    # Convert single item tuples to str
    for key, value in mapping.items():
        if len(value) == 1:
            mapping[key] = value[0]
    return mapping



def _recompute_node2coord(component: nx.DiGraph, quotient_relabeled: nx.DiGraph) -> dict:
    """Compute the new node2coord for the quotient graph"""
    logging.info("\tRecomputing node2coord")
    # Get the new_node2coord data
    old_node2coord = nx.get_node_attributes(G=component, name="coordinates")
    new_node2coord = {}
    for nodes in quotient_relabeled.nodes():
        if isinstance(nodes, str):  # node is untouched in quotient
            new_node2coord[nodes] = old_node2coord[nodes]
        elif isinstance(nodes, tuple):  # node was touched
            first_old_node = nodes[0]
            last_old_node = nodes[-1]
            transcript, start, _ = old_node2coord[first_old_node][0]
            _, _, end = old_node2coord[last_old_node][0]
            new_node2coord[nodes] = ((transcript, start, end),)
    return new_node2coord



def _recompute_edge2overlap(component: nx.DiGraph, quotient_relabeled: nx.DiGraph) -> dict:
    """Compute the new node2coord for the quotient graph"""
    logging.info("\tRecomputing edge2overlap")

    old_edge2overlap = nx.get_edge_attributes(G=component, name="overlaps")
    new_edge2overlap = {}

    for edge in quotient_relabeled.edges():
        node_u, node_v = edge

        if isinstance(node_u, tuple) and isinstance(node_v, tuple):
            new_edge2overlap[(node_u, node_v)] = old_edge2overlap[(node_u[-1], node_v[0])]

        elif isinstance(node_u, tuple) and isinstance(node_v, str):
            new_edge2overlap[(node_u, node_v)] = old_edge2overlap[(node_u[-1], node_v)]

        elif isinstance(node_u, str) and isinstance(node_v, tuple):
            new_edge2overlap[(node_u, node_v)] = old_edge2overlap[(node_u, node_v[0])]

        else:
            new_edge2overlap[(node_u, node_v)] = old_edge2overlap[(node_u, node_v)]

    return new_edge2overlap



def _compute_new_node_ids(quotient_relabeled: nx.DiGraph, component: nx.DiGraph) -> dict:
    """Compose the new node id for every collapsed node in nx.quotient_graph"""

    logging.info("\tRecomputing final node identifiers")

    quotient_mapping = {}
    old_node2coord = nx.get_node_attributes(G=component, name="coordinates")

    for node in quotient_relabeled.nodes():

        if isinstance(node, tuple):  # Collapsed node

            # Compute starting node and coordinates
            node = tuple(natsorted(node))
            first_node_id = node[0]
            last_node_id = node[-1]
            transcript_id, start, _ = old_node2coord[first_node_id][0]
            _, _, end = old_node2coord[last_node_id][0]

            # Compose new node
            new_node_id = "{transcript}:{start}-{end}".format(
                transcript=transcript_id,
                start=start,
                end=end
            )
            quotient_mapping[node] = new_node_id
        else:
            quotient_mapping[node] = node

    return quotient_mapping



def _sculpt_graph(splice_graph: nx.DiGraph, filled_edges: set) -> nx.DiGraph:
    """Apply sealer corrections in filled_edges to the splice graph"""

    logging.info("\tSculpting graph")

    if not filled_edges:
        return splice_graph

    # Compute the quotient graph
    def partition(node_u, node_v):
        """Function to test if node_u and node_v belong to the same partition of the graph"""
        graph = nx.DiGraph()
        graph.add_edges_from(filled_edges)
        if node_u in graph.nodes() and \
            node_v in graph.nodes() and \
            nx.has_path(G=graph, source=node_u, target=node_v):
            return True
        return False

    # Compute the quotient graph
    quotient = nx.quotient_graph(G=splice_graph, partition=partition)

    # Rename nodes (frozensets are tricky)
    mapping_sg_to_partition = _rename_nodes_from_collapse(quotient)
    quotient_relabeled = nx.relabel_nodes(
        G=quotient,
        mapping=mapping_sg_to_partition
    )

    # Recompute graph info
    node2coord = _recompute_node2coord(splice_graph, quotient_relabeled)
    edge2overlap = _recompute_edge2overlap(splice_graph, quotient_relabeled)
    node_ids = _compute_new_node_ids(component=splice_graph, quotient_relabeled=quotient_relabeled)

    # Set new info
    nx.set_node_attributes(G=quotient_relabeled, name="coordinates", values=node2coord)
    nx.set_edge_attributes(G=quotient_relabeled, name="overlaps", values=edge2overlap)
    component_final = nx.relabel_nodes(
        G=quotient_relabeled,
        mapping=node_ids
    )

    return component_final



def correct_splice_graph(splice_graph_dict: nx.DiGraph, args: dict) -> nx.DiGraph:
    """Try to correct small gaps and some overlaps (SNPs and indels) with abyss-sealer

    args = {
        "kmer": int,
        "max_gap_size": int,
        "input_bloom": str,
        "input_fasta": str
    }
    """
    logging.info("Correct splice graph with abyss-sealer")

    # Compose fasta with candidates to be filled
    sealer_input_fn = _prepare_sealer(
        splice_graph_dict=splice_graph_dict,
        args=args
    )

    # Run sealer
    sealer_output_fn = _run_sealer(
        sealer_input_fn=sealer_input_fn,
        args=args
    )

    # Collect results
    filled_edges = _collect_sealer_results(handle=sealer_output_fn)
    remove(sealer_input_fn)
    remove(sealer_output_fn)
    filled_edges_by_transcript = _filled_edges_by_transcript(filled_edges=filled_edges)
    del filled_edges

    # Complete the filled_edges_by_transcript dict
    for transcript in splice_graph_dict:
        if transcript not in filled_edges_by_transcript:
            filled_edges_by_transcript[transcript] = {}

    # Process each component separatedly
    for transcript_id, splice_graph in splice_graph_dict.items():
        logging.info("\tProcessing component %s", transcript_id)
        splice_graph_dict[transcript_id] = _sculpt_graph(
            splice_graph,
            filled_edges_by_transcript[transcript_id]
        )

    return splice_graph_dict
