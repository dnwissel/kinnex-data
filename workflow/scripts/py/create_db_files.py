import sys

import gffutils

with open(snakemake.log[0], "w") as f:
    sys.stderr = f
    sys.stdout = f

    # https://stackoverflow.com/questions/56402571/get-list-of-all-loaded-python-packages-versions-and-variables
    for module in sys.modules.values():
        if hasattr(module, "__version__"):
            print(module.__name__, module.__version__)
        else:
            print(module.__name__)

    db = gffutils.create.create_db(
        data=snakemake.input["input_gtf"],
        dbfn=snakemake.output["output_db"],
        force=True,
        verbose=False,
        force_gff=False,
        checklines=snakemake.params["checklines"],
        disable_infer_transcripts=True,
    )
