import re
import sys

import pandas as pd
import pysam

with open(snakemake.log[0], "w") as f:
    sys.stderr = f
    sys.stdout = f
    # https://stackoverflow.com/questions/56402571/get-list-of-all-loaded-python-packages-versions-and-variables
    for module in sys.modules.values():
        if hasattr(module, "__version__"):
            print(module.__name__, module.__version__)
        else:
            print(module.__name__)
    samfile = pysam.AlignmentFile(snakemake.input[0], "rb")

    # Adapted from: https://www.biostars.org/p/306041/
    read1 = None
    read2 = None
    read_name = []
    edit_pct = []
    splice_junctions = []
    for read in samfile:

        if (
            not read.is_paired
            or read.mate_is_unmapped
            or read.is_duplicate
            or read.is_secondary
            or read.is_supplementary
        ):
            continue

        if read.is_read2:
            read2 = read
        else:
            read1 = read
            read2 = None
            continue

        if (
            not read1 is None
            and not read2 is None
            and read1.query_name == read.query_name
        ):
            read_name += [read1.query_name]
            edit_pct += [
                (read1.get_tag("nM") + read2.get_tag("nM"))
                / (read1.query_alignment_length + read2.query_alignment_length)
            ]
            splice_junctions += [
                len(
                    list(
                        set(
                            re.findall(
                                snakemake.params["junction_regex"], read1.cigarstring
                            )
                            + re.findall(
                                snakemake.params["junction_regex"], read2.cigarstring
                            )
                        )
                    )
                )
            ]
    pd.DataFrame(
        {"qname": read_name, "edit_distance": edit_pct, "n_junctions": splice_junctions}
    ).to_csv(snakemake.output[0], sep="\t", index=False)
