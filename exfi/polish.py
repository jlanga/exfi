"""exfi.polish.py

exfi submodule to polish a bed4 dataframe by checking if in the overlap between
two exons there is the AG-GT splicing signal.
"""

import logging

def polish_bed4(bed4, transcriptome_dict):
    """
    Trim overlapping exons according to the AG-GT signal.
    """
    logging.info("Polishing BED4")

    polished = bed4.copy()

    logging.info("Get the transcript_id of the next exon")
    polished["chrom_next"] = polished["chrom"].shift(-1)

    logging.info('Get the name of the next exon')
    polished["name_next"] = polished["name"].shift(-1)

    logging.info('Get the start of the next exon')
    polished["chrom_start_next"] = polished["chrom_start"].shift(-1)

    logging.info('Get the end of the next exon')
    polished["chrom_end_next"] = polished["chrom_end"].shift(-1)

    logging.info('Remove rows with different transcripts')
    polished = polished\
        [polished["chrom"] == polished["chrom_next"]]

    logging.info('Cast from float to int (just in case)')
    polished = polished.astype({
        "chrom_start_next": int,
        "chrom_end_next": int
    })

    logging.info('Compute the overlap')
    polished["overlap"] = polished["chrom_end"] - polished["chrom_start_next"]

    logging.info('Throw away lines that cannot be polished')
    polished = polished[polished.overlap >= 4]

    logging.info('Get the entire transcript sequence')
    polished["sequence"] = polished.chrom.map(transcriptome_dict)

    logging.info('Prepare a column with the data required to extract the '
                 'overlapping sequence')
    polished["data_to_map"] = list(zip(
        polished.sequence,
        polished.chrom_start_next + 1,
        polished.chrom_end + 1
    ))

    logging.info('Get the overlapping sequence')
    polished["overlap_str"] = polished\
        .data_to_map\
        .map(lambda x: x[0][x[1]:x[2]])

    logging.info('Get the position in which the AGGT happens')
    polished["overlap_index"] = polished["overlap_str"].str.rfind("AGGT")

    logging.info('Throw away rows in which AGGT doesn\'t happen')
    polished = polished[polished.overlap_index >= 0]

    logging.info('Correct positions')
    polished["chrom_end_corrected"] = polished["chrom_end"] - 2
    polished["chrom_start_next_corrected"] = \
        polished["chrom_start_next"] + polished["overlap_index"] + 2

    logging.info('Organize the elements to correct')
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

    logging.info('Correct the starts')
    bed4_new.loc[starts_to_change.index.tolist()].chrom_start = \
        starts_to_change.chrom_start
    logging.info('Correct the ends')
    bed4_new.loc[ends_to_change.index.tolist()].chrom_end = \
        ends_to_change.chrom_end

    logging.info('Compose the new names')
    bed4_new = bed4_new.reset_index(drop=False)
    bed4_new["name"] = \
        bed4_new.chrom + ":" + \
        bed4_new.chrom_start.map(str) + "-" + \
        bed4_new.chrom_end.map(str)

    logging.info('Done')

    return bed4_new[["chrom", "chrom_start", "chrom_end", "name"]]
