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
    polished["chrom_next"] = polished["chrom"].shift(-1)

    # Get the name of the next exon
    polished["name_next"] = polished["name"].shift(-1)

    # Get the start of the next exon
    polished["chrom_start_next"] = polished["chrom_start"].shift(-1)

    # Get the end of the next exon
    polished["chrom_end_next"] = polished["chrom_end"].shift(-1)

    # Remove rows with different transcripts
    polished = polished\
        [polished["chrom"] == polished["chrom_next"]]

    # cast from float to int
    polished = polished.astype({"chrom_start_next": int, "chrom_end_next": int})

    # compute the overlap
    polished["overlap"] = polished["chrom_end"] - polished["chrom_start_next"]

    # Throw away lines that cannot be polished
    polished = polished[polished.overlap >= 4]

    # Get the entire transcript sequence
    polished["sequence"] = polished.chrom.map(transcriptome_dict)

    # Prepare a column with the data required to extract the overlapping seq
    polished["data_to_map"] = list(zip(
        polished.sequence,
        polished.chrom_start_next + 1,
        polished.chrom_end + 1
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
    polished["chrom_end_corrected"] = polished["chrom_end"] - 2
    polished["chrom_start_next_corrected"] = \
        polished["chrom_start_next"] + polished["overlap_index"] + 2

    # Organize the elements to correct
    ends_to_change = polished\
        [["name", "chrom_end_corrected"]]\
        .rename({"chrom_end_corrected": "chrom_end"}, axis=1)\
        .set_index("name")

    starts_to_change = polished\
        [["name_next", "chrom_start_next_corrected"]]\
        .rename(columns={
            "name_next": "name",
            "chrom_start_next_corrected":
            "chrom_start"
        })\
        .set_index("name")


    bed4_new = bed4.set_index("name")
    # Correct the starts
    bed4_new.loc[starts_to_change.index.tolist()].chrom_start = \
        starts_to_change.chrom_start
    # Correct the ends
    bed4_new.loc[ends_to_change.index.tolist()].chrom_end = \
        ends_to_change.chrom_end

    bed4_new = bed4_new.reset_index(drop=False)
    bed4_new["name"] = \
        bed4_new.chrom + ":" + \
        bed4_new.chrom_start.map(str) + "-" + \
        bed4_new.chrom_end.map(str)

    return bed4_new[["chrom", "chrom_start", "chrom_end", "name"]]
