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
    read_name = []
    tlen = []
    for ix, read in enumerate(samfile):
        if ix % 1000 == 0:
            print(ix)

        if read.is_duplicate or read.is_secondary or read.is_supplementary:
            continue
        read_name += [read.query_name]
        tlen += [read.template_length]
    pd.DataFrame({"qname": read_name, "tlen": tlen}).to_csv(
        snakemake.output[0], sep="\t", index=False
    )
