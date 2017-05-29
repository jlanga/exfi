#!/usr/bin/env python3

def _extract_loci_start_end(description):
    """(str) -> str, int, int
    Get coordinates from string
    """
    locus = description.split(":")[0]
    start, end = description.split(":")[1].split("-")
    return locus, int(start), int(end)

def _modify_description(description, bases_to_the_left, bases_to_the_right):
    """(str) -> str
    Process multiple descriptions at once
    """
    descriptions = description.split(" ")[1:] # First one is the exon_id
    modified_descriptions = []
    for description in descriptions:
        locus, start, end = _extract_loci_start_end(description)
        start += bases_to_the_left
        end -= bases_to_the_right
        modified_descriptions.append(
            "{locus}:{start}-{end}".format(
                locus = locus,
                start = start,
                end = end
            )
        )

    return " ".join(modified_descriptions)

def trim_sequence(iterable_records, bases_to_the_left=0, bases_to_the_right=0):
    """(iterable, int, int) -> generator
    Given a record, trim the bases to the left and the ones to the right.
    """

    iterable_records = (
        record for record in iterable_records \
        if len(record) > bases_to_the_left + bases_to_the_right
    )

    for record in iterable_records:
        # Trim bases
        n = len(record.seq)
        record.seq = record.seq[bases_to_the_left:n-bases_to_the_right]
        # Modify description
        record.description = _modify_description(
            record.description, bases_to_the_left, bases_to_the_right
        )
        yield record
