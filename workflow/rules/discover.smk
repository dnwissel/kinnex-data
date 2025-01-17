configfile: "config/config.yaml"


rule discover_run_stringtie2:
    input:
        reads="results/align/separate_sirv_and_gencode_illumina/{sample}/{sample}.aligned.gencode.sorted.bam",
        transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
    output:
        "results/discover/run_stringtie2/{sample}/transcriptome.gtf",
    log:
        "logs/discover/run_stringtie2/{sample}/log.out",
    conda:
        "../envs/standalone/stringtie2.yaml"
    threads: config["discover_isoforms_threads"]
    shell:
        """
        conda list > {log};
        stringtie -G {input.transcriptome} -p {threads} -o {output} {input.reads} &>> {log}
        """


rule discover_merge_stringtie2:
    input:
        discovered_transcriptomes=expand(
            "results/discover/run_stringtie2/{sample}/transcriptome.gtf",
            sample=config["sample_names"],
        ),
        reference_transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
    output:
        "results/discover/merge_stringtie2/transcriptome.gtf",
    log:
        "logs/discover/merge_stringtie2/log.out",
    conda:
        "../envs/standalone/stringtie2.yaml"
    shell:
        """
        conda list > {log};
        stringtie -G {input.reference_transcriptome} \
            -o {output} {input.discovered_transcriptomes} &>> {log}
        """


rule discover_run_bambu:
    input:
        reads=expand(
            "results/align/run_minimap2_gencode/{sample}/{sample}.aligned.sorted.bam",
            sample=config["sample_names"],
        ),
        transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    output:
        "results/discover/run_bambu/transcriptome.gtf",
    log:
        "logs/discover/run_bambu/log.out",
    conda:
        "../envs/r/bambu.yaml"
    threads: config["discover_isoforms_threads"]
    script:
        "../scripts/r/discover/run_bambu.R"


rule discover_fix_bambu_gene_ids:
    input:
        input_path="results/discover/run_bambu/transcriptome.gtf",
    output:
        output_path="results/discover/fix_bambu_gene_ids/transcriptome.gtf",
    log:
        "logs/discover/fix_bambu_gene_ids/log.out",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/bambu_fix_gene_ids.R"
