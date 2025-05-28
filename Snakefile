from snakemake.utils import min_version


configfile: "config/config.yaml"


min_version(config["snakemake_min_version"])


singularity: f"docker://condaforge/mambaforge:{config['mambaforge_version']}"


include: "workflow/rules/prepare.smk"
include: "workflow/rules/align.smk"
include: "workflow/rules/discover.smk"
include: "workflow/rules/downsample.smk"
include: "workflow/rules/quantify_subsampled.smk"
include: "workflow/rules/quantify.smk"
include: "workflow/rules/format_quantify.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/annotate.smk"
include: "workflow/rules/discover_sirv.smk"
include: "workflow/rules/plot_figures.smk"
include: "workflow/rules/make_tables.smk"


rule all:
    input:
        "results/plot_figures/figure_01.pdf",
        "results/plot_figures/figure_02.pdf",
        "results/plot_figures/figure_S01.pdf",
        "results/plot_figures/figure_S02.pdf",
        "results/plot_figures/figure_S03.pdf",
        "results/make_tables/table_S01.tsv",
        "results/make_tables/table_S02.tsv",
        "results/make_tables/table_S03.tsv",
        expand(
            "benchmarks/quantify_subsampled/{method}/{data_type}/{subsample_number}/{read_number}/{sample}.txt",
            method=[
                "run_salmon_illumina",
                "run_bambu",
                "run_oarfish",
                "run_kallisto_long",
                "run_isoquant",
                "run_salmon",
            ],
            data_type=["gencode"],
            subsample_number=["1"],
            read_number=["5000000.0", "10000000.0", "20000000.0", "30000000.0"],
            sample=[
                "day0-rep1",
                "day0-rep2",
                "day0-rep3",
                "day5-rep1",
                "day5-rep2",
                "day5-rep3",
            ],
        ),
        expand(
            "results/quantify_subsampled/run_oarfish_bootstrap/{data_type}/{subsample_number}/{read_number}/{sample}/{sample}.quant",
            data_type=["gencode"],
            subsample_number=["1"],
            read_number=["30000000.0"],
            sample=[
                "day0-rep1",
                "day0-rep2",
                "day0-rep3",
                "day5-rep1",
                "day5-rep2",
                "day5-rep3",
            ],
        ),
