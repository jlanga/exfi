#!/usr/bin/env python3

"""
exfi.io.read_gfa1: submodule to process a gfa1 into almost a splice graph
"""

import logging


def _overlap_str_to_int(overlap_str):
    """(str) -> int
    Modify overlap str to int:
    20G -> -20
    13M -> M
    """
    if not isinstance(overlap_str, str):
        raise TypeError("{overlap} is not str".format(overlap=overlap_str))
    letter = overlap_str[-1]
    if letter == "M":
        return int(overlap_str[:-1])
    elif letter == "G":
        return -int(overlap_str[:-1])
    else:
        raise ValueError("{letter} letter is not M or G".format(letter=letter))



def _process_segments(segments_raw):
    """(list of lists) -> dict

    Convert a list of ["S", node_id, str, *whatever] to a dict {node_id: str}
    """
    logging.info("\tProcessing segments")
    segments = {}
    for line in segments_raw:
        _, node_id, sequence, *_ = line
        segments[node_id] = sequence
    return segments



def _process_links(links_raw):
    """(list of lists) -> dict

    Convert a list of ["L", from, from_orient, to, to_ortient, overlap] to a dict
    {(from, to): overlap}
    """
    logging.info("\tProcessing links")
    links = {}
    for line in links_raw:
        _, node_u, _, node_v, _, overlap, *_ = line
        overlap = _overlap_str_to_int(overlap)
        links[(node_u, node_v)] = overlap
    return links



def _process_containments(containments_raw):
    """(list of lists) -> dict

    Convert a list of ["C", transcript_id, _, node_id, _, position, overlap] to a dict
    {node_id: (transcript_id, start, end)}
    """
    logging.info("\tProcessing containments")
    containments = {}
    for line in containments_raw:
        _, container, _, contained, _, position, overlap, *_ = line
        overlap = _overlap_str_to_int(overlap)
        start = int(position)
        end = start + overlap
        # containments[contained] = ((container, start, end), )
        if contained not in containments:
            containments[contained] = ()
        containments[contained] += ((container, start, end), )
    return containments



def _process_paths(containments_raw):
    """(list of lists) -> dict

    Conver a list of ["P", transcript_id, node1+,...,nodeN+]  to a dict {tra}
    """
    logging.info("\tProcessing paths")
    paths = {}
    for line in containments_raw:
        _, path_name, segment_names, *_ = line
        # Drop orientations!
        segment_names = tuple([segment_name[:-1] for segment_name in  segment_names.split(",")])
        paths[path_name] = segment_names
    return paths


def read_gfa1(filename):
    """(str) -> dict

    Convert GFA1 file to dict {
        "header": header list
        "segments": segment list,
        "links": link list,
        "cointainments": containments list,
        "paths": path list
    }
    """
    with open(filename, "r") as gfain:

        logging.info("Reading gfa1 {file}".format(file=filename))

        segments = []
        links = []
        containments = []
        paths = []

        for line in gfain:

            line = line.strip().split("\t")

            if line[0] == "H":
                header = line
            elif line[0] == "S":
                segments.append(line)
            elif line[0] == "L":
                links.append(line)
            elif line[0] == "C":
                containments.append(line)
            elif line[0] == "P":
                paths.append(line)

    return {
        "header": header[1],
        "segments": _process_segments(segments),
        "links": _process_links(links),
        "containments": _process_containments(containments),
        "paths": _process_paths(paths)
    }
