def _extract_loci_start_end(description):
    """(str) -> str, int, int
    Get coordinates from string
    """
    locus = description.split(":")[0]
    start, end = description.split(":")[1].split("-")
    return locus, int(start), int(end)


def _modify_description(
    description,
    bases_to_the_left=0,
    bases_to_the_right=0
):
    """(str) -> str
    Process multiple descriptions at once
    """
    descriptions = description.split(" ")[1:]  # First one is the exon_id
    modified_descriptions = []

    for description in descriptions:
        locus, start, end = _extract_loci_start_end(description)
        start += bases_to_the_left
        end -= bases_to_the_right
        modified_descriptions.append(
            "{locus}:{start}-{end}".format(
                locus=locus,
                start=start,
                end=end
            )
        )
    return " ".join(modified_descriptions)


def filter_by_extensibility(exons, bloom_filter, kmer):

    # Import shit
    from Bio import SeqIO  # To read and write fastas
    from exfi.filter_by_length import filter_by_length
    from exfi.extend import extend_left, extend_right
    from subprocess import Popen
    from itertools import chain
    import sys

    # Read raw exons. Use list since we will process it three times:
    # - one to compute left extensions,
    # - one to compute the right extensions, and
    # - one to trim and filter

    print("Computing all extensions", file=sys.stderr)
    # Compute left and right extensions
    # - Use `extend_left` and `extend_right`
    extensions = chain(
        extend_left(exons, kmer),
        extend_right(exons, kmer)
    )

    # - Write all of them in a temp file
    from tempfile import NamedTemporaryFile
    prefix = NamedTemporaryFile()
    extensions_raw = prefix.name + "_raw.fa"
    extensions_filtered = prefix.name + "_filtered.fa"
    SeqIO.write(
        sequences=extensions,
        handle=extensions_raw,
        format="fasta"
    )

    print("Checking which are possible", file=sys.stderr)
    # - Use `abyss-bloom kmers` to find out which extensions are there
    command_kmers = [  # Run abyss-bloom kmers
        "abyss-bloom", "kmers",
        "--kmer", str(kmer),
        "--verbose",
        bloom_filter,
        extensions_raw
    ]
    out_handle = open(extensions_filtered, "w")
    Popen(command_kmers, stdout=out_handle)
    # - Read the results (as a generator)
    extensions = SeqIO.parse(handle=extensions_filtered, format="fasta")
    print("Trimming by extensibility", file=sys.stderr)
    # - Process the extensions by name
    # i.e., from EXONXXXXXXXX[l|r]["ACGT"] just get the EXONXXXXXXX[l|r]
    extension_names = set(
        [extension.id.split(":")[0][:-1] for extension in extensions]
    )
    for exon in exons:
        if exon.id + "l" not in extension_names:
            exon.seq = exon.seq[1:]
            exon.description = _modify_description(
                exon.description,
                bases_to_the_left=1
            )
        if exon.id + "r" not in extension_names:
            exon.seq = exon.seq[:-1]
            exon.description = _modify_description(
                exon.description,
                bases_to_the_right=1
            )

    print("Filtering by length and writing to disk", file=sys.stderr)
    exons = filter_by_length(iterables=exons, length=kmer)

    for exon in exons:
        yield exon
