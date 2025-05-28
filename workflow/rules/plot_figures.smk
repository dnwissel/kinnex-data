configfile: "config/config.yaml"


singularity: f"docker://condaforge/mambaforge:{config['mambaforge_version']}"


rule plot_figures_install_ggpubfigs:
    output:
        "results/plot_figures/install_ggpubfigs/done.txt",
    params:
        version=config["ggpubfigs_version"],
    log:
        "logs/plot_figures/install_ggpubfigs/out.log",
    conda:
        "../envs/r/figures_tables.yaml"
    script:
        "../scripts/r/install_ggpubfigs.R"


rule plot_figures_make_tombstone:
    input:
        # PLotting extra packages
        "results/plot_figures/install_ggpubfigs/done.txt",
        # SQANTI
        "results/annotate/run_sqanti/kinnex_wtc_11/kinnex_wtc_11_classification.txt",
        # ORFanage
        "results/annotate/run_orfanage/transcriptome.gtf",
        # SIRV discovery
        expand(
            "results/discover_sirv/evaluate_performance/{data_type}_{rep}_{num_reads}",
            data_type=["pb", "illumina"],
            rep=["1"],
            num_reads=[
                "250000.0",
                "500000.0",
                "1000000.0",
                "2500000.0",
            ],
        ),
        # # SIRV quants
        expand(
            "results/format_quantify/{method}/sirv/transcript_counts_formatted.tsv",
            method=[
                "salmon",
                "salmon_illumina",
                "isoquant",
                "bambu",
                "oarfish",
                "kallisto",
            ],
        ),
        # GENCODE quants
        expand(
            "results/format_quantify/{method}/gencode/transcript_counts_formatted.tsv",
            method=[
                "salmon",
                "salmon_illumina",
                "isoquant",
                "bambu",
                "oarfish",
                "kallisto",
            ],
        ),
        # GENCODE + novel quantifications
        expand(
            "results/format_quantify/{method}_novel/gencode/transcript_counts_formatted.tsv",
            method=[
                "salmon",
                "salmon_illumina",
                "isoquant",
                "bambu",
                "oarfish",
                "kallisto",
            ],
        ),
        # Downsampled SIRV quantifications
        expand(
            "results/format_quantify_subsampled/{method}/sirv/{rep}/{depth}/transcript_counts_formatted.tsv",
            method=[
                "salmon",
                "salmon_illumina",
                "isoquant",
                "bambu",
                "oarfish",
                "kallisto",
            ],
            rep=[(i + 1) for i in range(5)],
            depth=["250000.0", "500000.0", "1000000.0", "2500000.0"],
        ),
        # Downsampled GENCODE quantifications
        expand(
            "results/format_quantify_subsampled/{method}/gencode/{rep}/{depth}/transcript_counts_formatted.tsv",
            method=[
                "salmon",
                "salmon_illumina",
                "isoquant",
                "bambu",
                "oarfish",
                "kallisto",
            ],
            rep=[1],
            depth=["5000000.0", "10000000.0", "20000000.0", "30000000.0"],
        ),
        # QC frames all
        expand(
            "results/qc/prepare_sampled_read_quality_frame_illumina/{sample}/{data_type}/quality_plot_frame.tsv",
            sample=[
                "day0-rep1",
                "day0-rep2",
                "day0-rep3",
                "day1-rep1",
                "day2-rep1",
                "day3-rep1",
                "day3-rep2",
                "day4-rep1",
                "day5-rep1",
                "day5-rep2",
                "day5-rep3",
            ],
            data_type=["gencode"],
        ),
        expand(
            "results/qc/prepare_sampled_read_quality_frame_kinnex/{sample}/{data_type}/quality_plot_frame.tsv",
            sample=[
                "day0-rep1",
                "day0-rep2",
                "day0-rep3",
                "day1-rep1",
                "day2-rep1",
                "day3-rep1",
                "day3-rep2",
                "day4-rep1",
                "day5-rep1",
                "day5-rep2",
                "day5-rep3",
            ],
            data_type=["gencode"],
        ),
        # Benchmarks
        expand(
            "benchmarks/quantify_subsampled/{method}/{data_type}/{subsample_number}/{read_number}/{sample}.txt",
            method=[
                "run_salmon_illumina",
                "run_bambu",
                "run_isoquant",
                "run_kallisto",
                "run_oarfish",
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
    output:
        "results/.figure_table_tombstone",
    shell:
        """
        touch {output}
        """


rule plot_figures_figure_01:
    input:
        "results/.figure_table_tombstone",
    output:
        "results/plot_figures/figure_01.pdf",
    log:
        "logs/plot_figures/figure_01/log.out",
    conda:
        "../envs/r/figures_tables.yaml"
    script:
        "../scripts/r/figure_01.R"


rule plot_figures_figure_02:
    input:
        "results/.figure_table_tombstone",
    output:
        "results/plot_figures/figure_02.pdf",
    params:
        depth="30000000.0",
    log:
        "logs/plot_figures/figure_02/log.out",
    conda:
        "../envs/r/figures_tables.yaml"
    script:
        "../scripts/r/figure_02.R"


rule plot_figures_figure_S01:
    input:
        "results/.figure_table_tombstone",
    output:
        "results/plot_figures/figure_S01.pdf",
    log:
        "logs/plot_figures/figure_S01/log.out",
    conda:
        "../envs/r/figures_tables.yaml"
    script:
        "../scripts/r/figure_S01.R"


rule plot_figures_figure_S02:
    input:
        "results/.figure_table_tombstone",
    output:
        "results/plot_figures/figure_S02.pdf",
    params:
        depth="30000000.0",
    log:
        "logs/plot_figures/figure_S02/log.out",
    conda:
        "../envs/r/figures_tables.yaml"
    script:
        "../scripts/r/figure_S02.R"


rule plot_figures_figure_S03:
    input:
        "results/.figure_table_tombstone",
    output:
        "results/plot_figures/figure_S03.pdf",
    params:
        depth="30000000.0",
    log:
        "logs/plot_figures/figure_S03/log.out",
    conda:
        "../envs/r/figures_tables.yaml"
    script:
        "../scripts/r/figure_S03.R"
