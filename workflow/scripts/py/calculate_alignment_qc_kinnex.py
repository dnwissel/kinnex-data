import re
import sys

import numpy as np
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
    read_name = []
    edit_pct = []
    splice_junctions = []
    for ix, read in enumerate(samfile):

        if read.is_duplicate or read.is_secondary or read.is_supplementary:
            continue
        read_name += [read.query_name]
        try:
            edit_pct += [(read.get_tag("NM") / read.query_alignment_length)]
        except KeyError:
            read_name = read_name[:-1]
            continue
        splice_junctions += [
            len(re.findall(snakemake.params["junction_regex"], read.cigarstring))
        ]
    pd.DataFrame(
        {"qname": read_name, "edit_distance": edit_pct, "n_junctions": splice_junctions}
    ).to_csv(snakemake.output[0], sep="\t", index=False)
