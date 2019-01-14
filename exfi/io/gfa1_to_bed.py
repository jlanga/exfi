#!/usr/bin/env python3

"""exfi.io.gfa1_to_bed.py: submodule to read a GFA1 file and convert it to BED4
"""

import pandas as pd
import numpy as np


def gfa1_to_bed4(filename):
    """Read a GFA1 file and convert it to BED4"""

    with open(filename, "r") as gfa:
        containments = pd.DataFrame(
            data=[
                x.strip().split("\t") for x in gfa.readlines() if x[0] == "C"
            ],
            columns=["RecordType", "Container", "ContainerOrient", "Contained",
                     "ContainedOrient", "Pos", "Overlap"],
            dtype=None
        )\
        .astype(dtype={"Pos": np.int})

    containments = containments.rename({
        "Container": "chrom",
        "Contained": "name"
    }, axis=1)
    containments["Overlap"] = containments["Overlap"]\
        .map(lambda x: int(x[:-1]))
    containments["chromStart"] = containments["Pos"]
    containments["chromEnd"] = containments["Pos"] + containments["Overlap"]
    containments = containments[["chrom", "chromStart", "chromEnd", "name"]]
    return containments
