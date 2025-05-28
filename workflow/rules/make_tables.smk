configfile: "config/config.yaml"


singularity: f"docker://condaforge/mambaforge:{config['mambaforge_version']}"


rule make_tables_table_S01:
    input:
        "results/.figure_table_tombstone",
    output:
        "results/make_tables/table_S01.tsv",
    log:
        "logs/make_tables/table_S01/log.out",
    conda:
        "../envs/r/figures_tables.yaml"
    script:
        "../scripts/r/table_S01.R"


rule make_tables_table_S02:
    input:
        "results/.figure_table_tombstone",
    output:
        "results/make_tables/table_S02.tsv",
    log:
        "logs/make_tables/table_S02/log.out",
    conda:
        "../envs/r/figures_tables.yaml"
    script:
        "../scripts/r/table_S02.R"


rule make_tables_table_S03:
    input:
        "results/.figure_table_tombstone",
    output:
        "results/make_tables/table_S03.tsv",
    log:
        "logs/make_tables/table_S03/log.out",
    conda:
        "../envs/r/figures_tables.yaml"
    script:
        "../scripts/r/table_S03.R"
