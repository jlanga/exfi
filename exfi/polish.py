"""exfi.polish.py

exfi submodule to polish a bed4 dataframe by checking if in the overlap between
two exons there is the AG-GT splicing signal.
"""

def polish_bed4(bed4, transcriptome_dict):
    """
    Trim overlapping exons according to the AG-GT signal.
    """
    polished = bed4.copy()

    # Get the transcript_id of the next exon
    polished["chromNext"] = polished["chrom"].shift(-1)

    # Get the name of the next exon
    polished["nameNext"] = polished["name"].shift(-1)

    # Get the start of the next exon
    polished["chromStartNext"] = polished["chromStart"].shift(-1)

    # Get the end of the next exon
    polished["chromEndNext"] = polished["chromEnd"].shift(-1)

    # Remove rows with different transcripts
    polished = polished\
        [polished["chrom"] == polished["chromNext"]]

    # cast from float to int
    polished = polished.astype({"chromStartNext": int, "chromEndNext": int})

    # compute the overlap
    polished["overlap"] = polished["chromEnd"] - polished["chromStartNext"]

    # Throw away lines that cannot be polished
    polished = polished[polished.overlap >= 4]

    # Get the entire transcript sequence
    polished["sequence"] = polished.chrom.map(transcriptome_dict)

    # Prepare a column with the data required to extract the overlapping seq
    polished["data_to_map"] = list(zip(
        polished.sequence,
        polished.chromStartNext + 1,
        polished.chromEnd + 1
    ))

    # Get the overlapping sequence
    polished["overlap_str"] = polished\
        .data_to_map\
        .map(lambda x: x[0][x[1]:x[2]])

    # Get the position in which the AGGT happens
    polished["overlap_index"] = polished["overlap_str"].str.rfind("AGGT")

    # Throw away rows in which AGGT doesn't happen
    polished = polished[polished.overlap_index >= 0]

    # Correct positions
    polished["chromEndCorrected"] = polished["chromEnd"] - 2
    polished["chromStartNextCorrected"] = \
        polished["chromStartNext"] + polished["overlap_index"] + 2

    # Organize the elements to correct
    ends_to_change = polished\
        [["name", "chromEndCorrected"]]\
        .rename({"chromEndCorrected": "chromEnd"}, axis=1)\
        .set_index("name")

    starts_to_change = polished\
        [["nameNext", "chromStartNextCorrected"]]\
        .rename(
            {"nameNext": "name", "chromStartNextCorrected": "chromStart"},
            axis=1
        )\
        .set_index("name")


    bed4_new = bed4.set_index("name")
    # Correct the starts
    bed4_new.loc[starts_to_change.index.tolist()].chromStart = \
        starts_to_change.chromStart
    # Correct the ends
    bed4_new.loc[ends_to_change.index.tolist()].chromEnd = \
        ends_to_change.chromEnd

    bed4_new = bed4_new.reset_index(drop=False)
    bed4_new["name"] = \
        bed4_new.chrom + ":" + \
        bed4_new.chromStart.map(str) + "-" + \
        bed4_new.chromEnd.map(str)

    return bed4_new[["chrom", "chromStart", "chromEnd", "name"]]
