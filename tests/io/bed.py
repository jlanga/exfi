#!/usr/bin/env python3

"""tests.io.bed.py: variables for testing bed dataframes"""

import pandas as pd
import numpy as np

from exfi.io.bed import \
    BED3_COLS, BED3_DTYPES, \
    BED4_COLS, BED4_DTYPES

BED3_EMPTY_FN = "tests/io/empty.bed"
BED3_SIMPLE_FN = "tests/io/simple.bed"
BED3_COMPLEX_FN = "tests/io/complex.bed"


BED3_EMPTY = pd.DataFrame(columns=BED3_COLS)
BED3_EMPTY = BED3_EMPTY.astype(BED3_DTYPES)

BED3_SIMPLE = pd.DataFrame(
    data=[("ENSDART00000161035.1", 0, 326)],
    columns=BED3_COLS
)

BED3_COMPLEX = pd.DataFrame(
    data=[
        ["ENSDART00000161035.1", 0, 326],
        ["ENSDART00000161035.1", 397, 472],
        ["ENSDART00000161035.1", 477, 523],
        ["ENSDART00000165342.1", 5, 127],
        ["ENSDART00000165342.1", 125, 304],
        ["ENSDART00000165342.1", 317, 460],
        ["ENSDART00000165342.1", 459, 592],
        ["ENSDART00000165342.1", 591, 650],
        ["ENSDART00000165342.1", 645, 746],
        ["ENSDART00000165342.1", 746, 851],
        ["ENSDART00000165342.1", 854, 886],
        ["ENSDART00000165342.1", 899, 953],
        ["ENSDART00000165342.1", 974, 1097],
        ["ENSDART00000165342.1", 1098, 1175],
        ["ENSDART00000165342.1", 1176, 1324]
    ],
    columns=BED3_COLS
)



BED4_EMPTY = pd.DataFrame(columns=BED4_COLS)
BED4_EMPTY = BED4_EMPTY.astype(BED4_DTYPES)


BED4_SIMPLE = pd.DataFrame(
    data=[("ENSDART00000161035.1", 0, 326, "ENSDART00000161035.1:0-326")],
    columns=BED4_COLS
)

BED4_COMPLEX = pd.DataFrame(
    data=[
        ["ENSDART00000161035.1", 0, 326, "ENSDART00000161035.1:0-326"],
        ["ENSDART00000161035.1", 397, 472, "ENSDART00000161035.1:397-472"],
        ["ENSDART00000161035.1", 477, 523, "ENSDART00000161035.1:477-523"],
        ["ENSDART00000165342.1", 5, 127, "ENSDART00000165342.1:5-127"],
        ["ENSDART00000165342.1", 125, 304, "ENSDART00000165342.1:125-304"],
        ["ENSDART00000165342.1", 317, 460, "ENSDART00000165342.1:317-460"],
        ["ENSDART00000165342.1", 459, 592, "ENSDART00000165342.1:459-592"],
        ["ENSDART00000165342.1", 591, 650, "ENSDART00000165342.1:591-650"],
        ["ENSDART00000165342.1", 645, 746, "ENSDART00000165342.1:645-746"],
        ["ENSDART00000165342.1", 746, 851, "ENSDART00000165342.1:746-851"],
        ["ENSDART00000165342.1", 854, 886, "ENSDART00000165342.1:854-886"],
        ["ENSDART00000165342.1", 899, 953, "ENSDART00000165342.1:899-953"],
        ["ENSDART00000165342.1", 974, 1097, "ENSDART00000165342.1:974-1097"],
        ["ENSDART00000165342.1", 1098, 1175, "ENSDART00000165342.1:1098-1175"],
        ["ENSDART00000165342.1", 1176, 1324, "ENSDART00000165342.1:1176-1324"]
    ],
    columns=BED4_COLS
)




NODE2COORDINATES_EMPTY = BED4_EMPTY.copy().set_index("name")

NODE2COORDINATES_SIMPLE = BED4_SIMPLE.copy().set_index("name")

NODE2COORDINATES_COMPLEX = BED4_COMPLEX.copy().set_index("name")



PATH2NODES_EMPTY = {}

PATH2NODES_SIMPLE = {
    "ENSDART00000161035.1" : ["ENSDART00000161035.1:0-326"]
}

PATH2NODES_COMPLEX = {
    "ENSDART00000161035.1": [
        "ENSDART00000161035.1:0-326", "ENSDART00000161035.1:397-472",
        "ENSDART00000161035.1:477-523"
    ],
    "ENSDART00000165342.1": [
        "ENSDART00000165342.1:5-127", "ENSDART00000165342.1:125-304",
        "ENSDART00000165342.1:317-460", "ENSDART00000165342.1:459-592",
        "ENSDART00000165342.1:591-650", "ENSDART00000165342.1:645-746",
        "ENSDART00000165342.1:746-851", "ENSDART00000165342.1:854-886",
        "ENSDART00000165342.1:899-953", "ENSDART00000165342.1:974-1097",
        "ENSDART00000165342.1:1098-1175", "ENSDART00000165342.1:1176-1324"
    ]
}



NODE2SEQUENCE_EMPTY = pd.DataFrame(columns=["name", "sequence"])

NODE2SEQUENCE_SIMPLE = pd.DataFrame(
    data=[
        ["ENSDART00000161035.1:0-326",
        "TGCACGGGTTTATTGTTCACAAAGAGATCGACAATGTGCGCAACTAAAATAAACATAGTACATTTTGATT"
        "ATACACGAACTTAAACTAAAGTCCAATCACACCTCCGCCCCGTTTCCACAGCAGCCTGTCAGGGTGGAGG"
        "AAAAGCGCGGCGGTCATGTGAGGCTCGAGCATCTCTCTCTCTCTCTCTCTCTCTCTCTCTACAGAATGAT"
        "AGAGGGAGCTCGTGAATCACATCATAGTCGTCCTCCCCTCATTCGTCCTCTCCAGCAGACACCGAAAAAC"
        "TGCGTTCATGCCAAAATGGGATGTGGAAATTCCTCCGCCACGAGCA"]
    ],
    columns=["name", "sequence"]
)

NODE2SEQUENCE_COMPLEX = pd.DataFrame(
    data=[
        [
            "ENSDART00000161035.1:0-326",
            "TGCACGGGTTTATTGTTCACAAAGAGATCGACAATGTGCGCAACTAAAATAAACATAGTACATTTT"
            "GATTATACACGAACTTAAACTAAAGTCCAATCACACCTCCGCCCCGTTTCCACAGCAGCCTGTCAG"
            "GGTGGAGGAAAAGCGCGGCGGTCATGTGAGGCTCGAGCATCTCTCTCTCTCTCTCTCTCTCTCTCT"
            "CTACAGAATGATAGAGGGAGCTCGTGAATCACATCATAGTCGTCCTCCCCTCATTCGTCCTCTCCA"
            "GCAGACACCGAAAAACTGCGTTCATGCCAAAATGGGATGTGGAAATTCCTCCGCCACGAGCA"
        ], [
            "ENSDART00000161035.1:397-472",
            "AGGAACTACGGTGGAGTGTATGTGGGTCTTCCTGCTGATCTGACTGCAGTCGCTGCCAGTCAGTCC"
            "AAATCAACA"
        ], [
            "ENSDART00000161035.1:477-523",
            "AGTCAACAGATGTTTATTGCAGACCTTCAGATAAAACAACATAGAA"
        ], [
            "ENSDART00000165342.1:5-127",
            "TGGAGCTGAAGCCGAGTATCTTGGTATTGGACTGGAACAGAAATCCAGCAAAAACTTTAAGGGAAA"
            "TCACTTTCATTTCATGATCGAAAAACTCCCGCAGATCATAAAAGAGTGGAAGGAAG"
        ], [
            "ENSDART00000165342.1:125-304",
            "AGGACCTGTAGTAGAAACAAAACTAGGATCTCTGAGAGGTGCCTTCTTGACTGTGAAGGGCAAGGA"
            "CACAATAGTCAATAGTTATCTAGGTGTGCCGTTCGCCAAGCCGCCTGTAGGACCCCTGAGACTTGC"
            "TCGACCACAGGCTGCAGAGAAATGGCAAGGAGTTAGAGATGCCACCA"
        ], [
            "ENSDART00000165342.1:317-460",
            "GTGCCTCCAGGAAAGGCAAATGACTGTAACTGAACTGGAGTTTCTATCGATGGATGTGGAGGTTCC"
            "TGAGGTCTCGGAGGATTGCCTGTATCTTAACATCTACACCCCAGTTAAACCTGGACAAGGAGACAA"
            "GAAGTTACCAG"
        ], [
            "ENSDART00000165342.1:459-592",
            "GTCATGGTTTGGATTCATGGTGGAGGACTCTCTCTTGGATCGGCTTCAATGTATGATGGCTCTGTT"
            "CTGGCTGCGTATCAGGATGTGGTCGTGGTGCTCATTCAGTACAGATTGGGTCTTCTGGGGTTCTTA"
            "A"
        ], [
            "ENSDART00000165342.1:591-650",
            "AGCACCGGAGACGAGCATGCGCCAGGAAACTATGGTTTTCTGGATCAAGTAGCTGCCCT"
        ], [
            "ENSDART00000165342.1:645-746",
            "GCCCTTCAGTGGGTTCAGGAGAACATCCACAGCTTCGGTGGAGATCCTGGATCAGTGACCATCTTT"
            "GGAGAGTCTGCTGGAGGAATCAGTGTATCCACGCT"
        ], [
            "ENSDART00000165342.1:746-851",
            "GATTCTTTCCCCGCTGGCGTCTGGACTGTTTCATCGCGCCATTGCAGAAAGTGGAACTGCCTTCTG"
            "GGATGGTTTAGTCATGGCTGATCCTTTTCAGAGAGCCCA"
        ], [
            "ENSDART00000165342.1:854-886",
            "TGCAGCCAAACAATGCAACTGTGACAGCAGCA",
        ], [
            "ENSDART00000165342.1:899-953",
            "TGTCGACTGCATTATGCACTGGTCTGAAGAGGAGGCTCTGGAATGTGCTAAAAA"
        ], [
            "ENSDART00000165342.1:974-1097",
            "CGTTGCTGTAGATTCTTATTTCCTTCCCAAACCCATCGAGGAGATTGTTGAGAAACAAGAGTTTAG"
            "TAAAGTTCCTCTCATCAACGGCATTAACAATGATGAGTTTGGCTTCTTGTTGGCTGA"
        ], [
            "ENSDART00000165342.1:1098-1175",
            "TATTTCTTGGGTCCTGAATGGATGAATGGGTTGAAAAGAGAGCAAATCGCTGAAGCCTTGACGCTC"
            "ACATATCCTGA"
        ], [
            "ENSDART00000165342.1:1176-1324",
            "CCCAAGGATCGATGGATCATTGATCTGGTGGCGAAGGAATATCTGGGCGACACACACGACCCCATT"
            "GAAATCCGTGAAGTTTATCGGGAGATGATGGGAGACGTGCTGTTTAACATCCCTGCCCTGCAACTG"
            "GCAAAACACCACAGCG"
        ]
    ],
    columns=["name", "sequence"]
)



EDGE2OVERLAP_EMPTY = pd.DataFrame(columns=["u", "v", "overlap"])
EDGE2OVERLAP_EMPTY = EDGE2OVERLAP_EMPTY.astype({"overlap": np.int64})

EDGE2OVERLAP_SIMPLE = pd.DataFrame(columns=["u", "v", "overlap"])
EDGE2OVERLAP_SIMPLE = EDGE2OVERLAP_SIMPLE.astype({"overlap": np.int64})

EDGE2OVERLAP_COMPLEX = pd.DataFrame(
    data=[
        ["ENSDART00000161035.1:0-326", "ENSDART00000161035.1:397-472", -71],
        ["ENSDART00000161035.1:397-472", "ENSDART00000161035.1:477-523", -5],
        ["ENSDART00000165342.1:5-127", "ENSDART00000165342.1:125-304", 2],
        ["ENSDART00000165342.1:125-304", "ENSDART00000165342.1:317-460", -13],
        ["ENSDART00000165342.1:317-460", "ENSDART00000165342.1:459-592", 1],
        ["ENSDART00000165342.1:459-592", "ENSDART00000165342.1:591-650", 1],
        ["ENSDART00000165342.1:591-650", "ENSDART00000165342.1:645-746", 5],
        ["ENSDART00000165342.1:645-746", "ENSDART00000165342.1:746-851", 0],
        ["ENSDART00000165342.1:746-851", "ENSDART00000165342.1:854-886", -3],
        ["ENSDART00000165342.1:854-886", "ENSDART00000165342.1:899-953", -13],
        ["ENSDART00000165342.1:899-953", "ENSDART00000165342.1:974-1097", -21],
        ["ENSDART00000165342.1:974-1097", "ENSDART00000165342.1:1098-1175", -1],
        ["ENSDART00000165342.1:1098-1175", "ENSDART00000165342.1:1176-1324", -1]
    ],
    columns=["u", "v", "overlap"]
)
EDGE2OVERLAP_COMPLEX = EDGE2OVERLAP_COMPLEX.astype({"overlap": np.int64})



BED3_ENSEMBL = pd.DataFrame(
    data=[
        ['ENSDART00000161842', 0, 43],
        ['ENSDART00000161842', 43, 192],
        ['ENSDART00000161842', 192, 351],
        ['ENSDART00000161842', 351, 417],
        ['ENSDART00000161842', 417, 531],
        ['ENSDART00000161842', 531, 695],
        ['ENSDART00000161842', 695, 722],
        ['ENSDART00000165461', 0, 236],
        ['ENSDART00000165461', 236, 347],
        ['ENSDART00000166393', 0, 213],
        ['ENSDART00000166393', 213, 362],
        ['ENSDART00000166393', 362, 521],
        ['ENSDART00000166393', 521, 623],
        ['ENSDART00000166393', 623, 737],
        ['ENSDART00000166393', 737, 901],
        ['ENSDART00000166393', 901, 1060],
        ['ENSDART00000166393', 1060, 1173],
        ['ENSDART00000166393', 1173, 1612],
        ['ENSDART00000170165', 0, 68],
        ['ENSDART00000170165', 68, 217],
        ['ENSDART00000170165', 217, 398],
        ['ENSDART00000170165', 398, 595],
        ['ENSDART00000170165', 595, 759],
        ['ENSDART00000170165', 759, 918],
        ['ENSDART00000170165', 918, 1031],
        ['ENSDART00000170165', 1031, 1470],
        ['ENSDART00000170877', 0, 141],
        ['ENSDART00000170877', 141, 290],
        ['ENSDART00000170877', 290, 449],
        ['ENSDART00000170877', 449, 551],
        ['ENSDART00000170877', 551, 590],
        ['ENSDART00000171631', 0, 90],
        ['ENSDART00000157701', 0, 362],
        ['ENSDART00000157701', 362, 426],
        ['ENSDART00000157701', 426, 545],
        ['ENSDART00000158290', 0, 444],
        ['ENSDART00000158290', 444, 563],
        ['ENSDART00000164359', 0, 277],
        ['ENSDART00000164359', 277, 353],
        ['ENSDART00000164359', 353, 464],
        ['ENSDART00000164359', 464, 601],
        ['ENSDART00000164359', 601, 665],
        ['ENSDART00000164359', 665, 1018],
        ['ENSDART00000167898', 0, 176],
        ['ENSDART00000167898', 176, 287],
        ['ENSDART00000167898', 287, 424],
        ['ENSDART00000167898', 424, 488],
        ['ENSDART00000167898', 488, 605]
    ],
    columns=BED3_COLS
)
BED3_GMAP = pd.DataFrame(
    data=[
        ['ENSDART00000171570', 0, 61],
        ['ENSDART00000171570', 61, 284],
        ['ENSDART00000171570', 284, 462],
        ['ENSDART00000171570', 462, 571],
        ['ENSDART00000171570', 571, 673],
        ['ENSDART00000171570', 673, 766],
        ['ENSDART00000171570', 766, 934],
        ['ENSDART00000171570', 934, 1024],
        ['ENSDART00000171570', 1024, 1133],
        ['ENSDART00000171570', 1133, 1289],
        ['ENSDART00000171570', 1289, 1371],
        ['ENSDART00000157830', 0, 90],
        ['ENSDART00000157830', 90, 431],
        ['ENSDART00000158772', 0, 49],
        ['ENSDART00000158772', 49, 348],
        ['ENSDART00000158814', 0, 64],
        ['ENSDART00000158814', 64, 396],
        ['ENSDART00000159795', 0, 46],
        ['ENSDART00000159795', 46, 419],
        ['ENSDART00000160202', 0, 30],
        ['ENSDART00000160202', 30, 378],
        ['ENSDART00000160762', 0, 43],
        ['ENSDART00000160762', 43, 345],
        ['ENSDART00000160996', 0, 40],
        ['ENSDART00000160996', 40, 369],
        ['ENSDART00000161368', 0, 43],
        ['ENSDART00000161368', 43, 354],
        ['ENSDART00000161463', 0, 43],
        ['ENSDART00000161463', 43, 417],
        ['ENSDART00000162456', 0, 43],
        ['ENSDART00000162456', 43, 366],
        ['ENSDART00000163675', 0, 43],
        ['ENSDART00000163675', 43, 339],
        ['ENSDART00000163851', 0, 40],
        ['ENSDART00000163851', 40, 330],
        ['ENSDART00000164110', 0, 67],
        ['ENSDART00000164110', 67, 399],
        ['ENSDART00000164309', 0, 46],
        ['ENSDART00000164309', 46, 390],
        ['ENSDART00000164489', 0, 40],
        ['ENSDART00000164489', 40, 350],
        ['ENSDART00000164491', 0, 43],
        ['ENSDART00000164491', 43, 366],
        ['ENSDART00000165410', 0, 52],
        ['ENSDART00000165410', 52, 350],
        ['ENSDART00000166029', 0, 61],
        ['ENSDART00000166029', 61, 381],
        ['ENSDART00000166882', 0, 43],
        ['ENSDART00000166882', 43, 420],
        ['ENSDART00000166892', 0, 332],
        ['ENSDART00000167404', 0, 40],
        ['ENSDART00000167404', 40, 357],
        ['ENSDART00000167409', 0, 43],
        ['ENSDART00000167409', 43, 420],
        ['ENSDART00000167805', 0, 40],
        ['ENSDART00000167805', 40, 378],
        ['ENSDART00000168039', 0, 90],
        ['ENSDART00000168039', 90, 470],
        ['ENSDART00000170399', 0, 70],
        ['ENSDART00000170399', 70, 380],
        ['ENSDART00000170523', 0, 55],
        ['ENSDART00000170804', 0, 43],
        ['ENSDART00000170804', 43, 366],
        ['ENSDART00000171020', 0, 40],
        ['ENSDART00000171020', 40, 383],
        ['ENSDART00000171201', 0, 40],
        ['ENSDART00000171201', 40, 397],
        ['ENSDART00000171344', 0, 106],
        ['ENSDART00000171344', 106, 462],
        ['ENSDART00000171772', 0, 90],
        ['ENSDART00000171772', 90, 483],
        ['ENSDART00000172037', 0, 43],
        ['ENSDART00000172037', 43, 344],
        ['ENSDART00000172182', 0, 61],
        ['ENSDART00000172182', 61, 413],
        ['ENSDART00000172374', 0, 355]
    ],
    columns=BED3_COLS
)



BED4_SIMPLE_POLISHED = BED4_SIMPLE

BED4_COMPLEX_POLISHED = pd.DataFrame(
    data=[
        ["ENSDART00000161035.1", 0, 326, "ENSDART00000161035.1:0-326"],
        ["ENSDART00000161035.1", 397, 472, "ENSDART00000161035.1:397-472"],
        ["ENSDART00000161035.1", 477, 523, "ENSDART00000161035.1:477-523"],
        ["ENSDART00000165342.1", 5, 127, "ENSDART00000165342.1:5-127"],
        ["ENSDART00000165342.1", 125, 304, "ENSDART00000165342.1:125-304"],
        ["ENSDART00000165342.1", 317, 460, "ENSDART00000165342.1:317-460"],
        ["ENSDART00000165342.1", 459, 592, "ENSDART00000165342.1:459-592"],
        ["ENSDART00000165342.1", 591, 650, "ENSDART00000165342.1:591-650"],
        ["ENSDART00000165342.1", 645, 746, "ENSDART00000165342.1:645-746"],
        ["ENSDART00000165342.1", 746, 851, "ENSDART00000165342.1:746-851"],
        ["ENSDART00000165342.1", 854, 886, "ENSDART00000165342.1:854-886"],
        ["ENSDART00000165342.1", 899, 953, "ENSDART00000165342.1:899-953"],
        ["ENSDART00000165342.1", 974, 1097, "ENSDART00000165342.1:974-1097"],
        ["ENSDART00000165342.1", 1098, 1175, "ENSDART00000165342.1:1098-1175"],
        ["ENSDART00000165342.1", 1176, 1324, "ENSDART00000165342.1:1176-1324"]
    ],
    columns=BED4_COLS
)



NODE2SEQUENCE_COMPLEX_SOFT = pd.DataFrame(
    data=[
        [
            "ENSDART00000161035.1:0-326",
            "TGCACGGGTTTATTGTTCACAAAGAGATCGACAATGTGCGCAACTAAAATAAACATAGTACATTTT"
            "GATTATACACGAACTTAAACTAAAGTCCAATCACACCTCCGCCCCGTTTCCACAGCAGCCTGTCAG"
            "GGTGGAGGAAAAGCGCGGCGGTCATGTGAGGCTCGAGCATCTCTCTCTCTCTCTCTCTCTCTCTCT"
            "CTACAGAATGATAGAGGGAGCTCGTGAATCACATCATAGTCGTCCTCCCCTCATTCGTCCTCTCCA"
            "GCAGACACCGAAAAACTGCGTTCATGCCAAAATGGGATGTGGAAATTCCTCCGCCACGAGCA"
        ], [
            "ENSDART00000161035.1:397-472",
            "AGGAACTACGGTGGAGTGTATGTGGGTCTTCCTGCTGATCTGACTGCAGTCGCTGCCAGTCAGTCC"
            "AAATCAACA"
        ], [
            "ENSDART00000161035.1:477-523",
            "AGTCAACAGATGTTTATTGCAGACCTTCAGATAAAACAACATAGAA"
        ], [
            "ENSDART00000165342.1:5-127",
            "TGGAGCTGAAGCCGAGTATCTTGGTATTGGACTGGAACAGAAATCCAGCAAAAACTTTAAGGGAAA"
            "TCACTTTCATTTCATGATCGAAAAACTCCCGCAGATCATAAAAGAGTGGAAGGAag"
        ], [
            "ENSDART00000165342.1:125-304",
            "agGACCTGTAGTAGAAACAAAACTAGGATCTCTGAGAGGTGCCTTCTTGACTGTGAAGGGCAAGGA"
            "CACAATAGTCAATAGTTATCTAGGTGTGCCGTTCGCCAAGCCGCCTGTAGGACCCCTGAGACTTGC"
            "TCGACCACAGGCTGCAGAGAAATGGCAAGGAGTTAGAGATGCCACCA"
        ], [
            "ENSDART00000165342.1:317-460",
            "GTGCCTCCAGGAAAGGCAAATGACTGTAACTGAACTGGAGTTTCTATCGATGGATGTGGAGGTTCC"
            "TGAGGTCTCGGAGGATTGCCTGTATCTTAACATCTACACCCCAGTTAAACCTGGACAAGGAGACAA"
            "GAAGTTACCAg"
        ], [
            "ENSDART00000165342.1:459-592",
            "gTCATGGTTTGGATTCATGGTGGAGGACTCTCTCTTGGATCGGCTTCAATGTATGATGGCTCTGTT"
            "CTGGCTGCGTATCAGGATGTGGTCGTGGTGCTCATTCAGTACAGATTGGGTCTTCTGGGGTTCTTA"
            "a"
        ], [
            "ENSDART00000165342.1:591-650",
            "aGCACCGGAGACGAGCATGCGCCAGGAAACTATGGTTTTCTGGATCAAGTAGCTgccct"
        ], [
            "ENSDART00000165342.1:645-746",
            "gccctTCAGTGGGTTCAGGAGAACATCCACAGCTTCGGTGGAGATCCTGGATCAGTGACCATCTTT"
            "GGAGAGTCTGCTGGAGGAATCAGTGTATCCACGCT"
        ], [
            "ENSDART00000165342.1:746-851",
            "GATTCTTTCCCCGCTGGCGTCTGGACTGTTTCATCGCGCCATTGCAGAAAGTGGAACTGCCTTCTG"
            "GGATGGTTTAGTCATGGCTGATCCTTTTCAGAGAGCCCA"
        ], [
            "ENSDART00000165342.1:854-886",
            "TGCAGCCAAACAATGCAACTGTGACAGCAGCA",
        ], [
            "ENSDART00000165342.1:899-953",
            "TGTCGACTGCATTATGCACTGGTCTGAAGAGGAGGCTCTGGAATGTGCTAAAAA"
        ], [
            "ENSDART00000165342.1:974-1097",
            "CGTTGCTGTAGATTCTTATTTCCTTCCCAAACCCATCGAGGAGATTGTTGAGAAACAAGAGTTTAG"
            "TAAAGTTCCTCTCATCAACGGCATTAACAATGATGAGTTTGGCTTCTTGTTGGCTGA"
        ], [
            "ENSDART00000165342.1:1098-1175",
            "TATTTCTTGGGTCCTGAATGGATGAATGGGTTGAAAAGAGAGCAAATCGCTGAAGCCTTGACGCTC"
            "ACATATCCTGA"
        ], [
            "ENSDART00000165342.1:1176-1324",
            "CCCAAGGATCGATGGATCATTGATCTGGTGGCGAAGGAATATCTGGGCGACACACACGACCCCATT"
            "GAAATCCGTGAAGTTTATCGGGAGATGATGGGAGACGTGCTGTTTAACATCCCTGCCCTGCAACTG"
            "GCAAAACACCACAGCG"
        ]
    ],
    columns=["name", "sequence"]
)

NODE2SEQUENCE_COMPLEX_HARD = pd.DataFrame(
    data=[
        [
            "ENSDART00000161035.1:0-326",
            "TGCACGGGTTTATTGTTCACAAAGAGATCGACAATGTGCGCAACTAAAATAAACATAGTACATTTT"
            "GATTATACACGAACTTAAACTAAAGTCCAATCACACCTCCGCCCCGTTTCCACAGCAGCCTGTCAG"
            "GGTGGAGGAAAAGCGCGGCGGTCATGTGAGGCTCGAGCATCTCTCTCTCTCTCTCTCTCTCTCTCT"
            "CTACAGAATGATAGAGGGAGCTCGTGAATCACATCATAGTCGTCCTCCCCTCATTCGTCCTCTCCA"
            "GCAGACACCGAAAAACTGCGTTCATGCCAAAATGGGATGTGGAAATTCCTCCGCCACGAGCA"
        ], [
            "ENSDART00000161035.1:397-472",
            "AGGAACTACGGTGGAGTGTATGTGGGTCTTCCTGCTGATCTGACTGCAGTCGCTGCCAGTCAGTCC"
            "AAATCAACA"
        ], [
            "ENSDART00000161035.1:477-523",
            "AGTCAACAGATGTTTATTGCAGACCTTCAGATAAAACAACATAGAA"
        ], [
            "ENSDART00000165342.1:5-127",
            "TGGAGCTGAAGCCGAGTATCTTGGTATTGGACTGGAACAGAAATCCAGCAAAAACTTTAAGGGAAA"
            "TCACTTTCATTTCATGATCGAAAAACTCCCGCAGATCATAAAAGAGTGGAAGGANN"
        ], [
            "ENSDART00000165342.1:125-304",
            "NNGACCTGTAGTAGAAACAAAACTAGGATCTCTGAGAGGTGCCTTCTTGACTGTGAAGGGCAAGGA"
            "CACAATAGTCAATAGTTATCTAGGTGTGCCGTTCGCCAAGCCGCCTGTAGGACCCCTGAGACTTGC"
            "TCGACCACAGGCTGCAGAGAAATGGCAAGGAGTTAGAGATGCCACCA"
        ], [
            "ENSDART00000165342.1:317-460",
            "GTGCCTCCAGGAAAGGCAAATGACTGTAACTGAACTGGAGTTTCTATCGATGGATGTGGAGGTTCC"
            "TGAGGTCTCGGAGGATTGCCTGTATCTTAACATCTACACCCCAGTTAAACCTGGACAAGGAGACAA"
            "GAAGTTACCAN"
        ], [
            "ENSDART00000165342.1:459-592",
            "NTCATGGTTTGGATTCATGGTGGAGGACTCTCTCTTGGATCGGCTTCAATGTATGATGGCTCTGTT"
            "CTGGCTGCGTATCAGGATGTGGTCGTGGTGCTCATTCAGTACAGATTGGGTCTTCTGGGGTTCTTA"
            "N"
        ], [
            "ENSDART00000165342.1:591-650",
            "NGCACCGGAGACGAGCATGCGCCAGGAAACTATGGTTTTCTGGATCAAGTAGCTNNNNN"
        ], [
            "ENSDART00000165342.1:645-746",
            "NNNNNTCAGTGGGTTCAGGAGAACATCCACAGCTTCGGTGGAGATCCTGGATCAGTGACCATCTTT"
            "GGAGAGTCTGCTGGAGGAATCAGTGTATCCACGCT"
        ], [
            "ENSDART00000165342.1:746-851",
            "GATTCTTTCCCCGCTGGCGTCTGGACTGTTTCATCGCGCCATTGCAGAAAGTGGAACTGCCTTCTG"
            "GGATGGTTTAGTCATGGCTGATCCTTTTCAGAGAGCCCA"
        ], [
            "ENSDART00000165342.1:854-886",
            "TGCAGCCAAACAATGCAACTGTGACAGCAGCA",
        ], [
            "ENSDART00000165342.1:899-953",
            "TGTCGACTGCATTATGCACTGGTCTGAAGAGGAGGCTCTGGAATGTGCTAAAAA"
        ], [
            "ENSDART00000165342.1:974-1097",
            "CGTTGCTGTAGATTCTTATTTCCTTCCCAAACCCATCGAGGAGATTGTTGAGAAACAAGAGTTTAG"
            "TAAAGTTCCTCTCATCAACGGCATTAACAATGATGAGTTTGGCTTCTTGTTGGCTGA"
        ], [
            "ENSDART00000165342.1:1098-1175",
            "TATTTCTTGGGTCCTGAATGGATGAATGGGTTGAAAAGAGAGCAAATCGCTGAAGCCTTGACGCTC"
            "ACATATCCTGA"
        ], [
            "ENSDART00000165342.1:1176-1324",
            "CCCAAGGATCGATGGATCATTGATCTGGTGGCGAAGGAATATCTGGGCGACACACACGACCCCATT"
            "GAAATCCGTGAAGTTTATCGGGAGATGATGGGAGACGTGCTGTTTAACATCCCTGCCCTGCAACTG"
            "GCAAAACACCACAGCG"
        ]
    ],
    columns=["name", "sequence"]
)