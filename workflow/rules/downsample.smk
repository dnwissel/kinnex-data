configfile: "config/config.yaml"


singularity: f"docker://condaforge/mambaforge:{config['mambaforge_version']}"


rule downsample_long_read_fastq:
    input:
        "results/prepare/convert_mapped_bams_to_fastq/{sample}/{type}.reads.fastq.gz",
    output:
        "results/downsample/convert_mapped_bams_to_fastq_long_reads/{subsample_number}_{read_number}_{type}/{sample}.fastq.gz",
    params:
        seed=lambda wc: wc.subsample_number,
        number_to_sample=lambda wc: wc.read_number,
    threads: config["stall_io_threads"]
    log:
        "logs/downsample/convert_mapped_bams_to_fastq/sirv/{subsample_number}/{read_number}/{type}/{sample}.log",
    conda:
        "../envs/standalone/seqtk.yaml"
    shell:
        """
        conda list &> {log};
        seqtk sample -2 -s {params.seed} {input} \
            {params.number_to_sample} | gzip > {output} 2>> {log}
        """


rule downsample_minimap2_transcriptome:
    input:
        input_name="results/align/run_minimap2_transcriptome_{type}/{sample}/{sample}.aligned.bam",
    output:
        output_name="results/downsample/run_minimap2_transcriptome_{type}/{subsample_number}/{read_number}/{sample}/{sample}.aligned.bam",
    params:
        seed=lambda wc: wc.subsample_number,
        number_to_sample=lambda wc: wc.read_number,
        is_transcriptome=True,
        read_name_sam=config["read_name_sam"],
    threads: config["stall_io_threads"]
    log:
        "logs/downsample/run_minimap2_transcriptome_{type}/{subsample_number}/{read_number}/{sample}.log",
    conda:
        "../envs/r/rsamtools.yaml"
    script:
        "../scripts/r/subsample_bam.R"


rule downsample_minimap2:
    input:
        input_name="results/align/run_minimap2_{type}/{sample}/{sample}.aligned.sorted.bam",
    output:
        output_name="results/downsample/run_minimap2_{type}/{subsample_number}/{read_number}/{sample}/{sample}.aligned.sorted.bam",
    params:
        seed=lambda wc: wc.subsample_number,
        number_to_sample=lambda wc: wc.read_number,
        is_transcriptome=False,
        read_name_sam=config["read_name_sam"],
    threads: config["stall_io_threads"]
    log:
        "logs/downsample/run_minimap2_{type}/{subsample_number}/{read_number}/{sample}.log",
    conda:
        "../envs/r/rsamtools.yaml"
    script:
        "../scripts/r/subsample_bam.R"


rule downsample_short_read_fastq:
    input:
        first_reads="results/align/convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.reads_1.fastq.gz",
        second_reads="results/align/convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.reads_2.fastq.gz",
    output:
        first_reads="results/downsample/convert_mapped_bams_to_fastq/{subsample_number}_{read_number}_{type}/{sample}-r1.fastq.gz",
        second_reads="results/downsample/convert_mapped_bams_to_fastq/{subsample_number}_{read_number}_{type}/{sample}-r2.fastq.gz",
    params:
        seed=lambda wc: wc.subsample_number,
        number_to_sample=lambda wc: wc.read_number,
    threads: config["stall_io_threads"]
    log:
        "logs/downsample/convert_mapped_bams_to_fastq/{type}/{subsample_number}/{read_number}/{sample}.log",
    conda:
        "../envs/standalone/seqtk.yaml"
    shell:
        """
        conda list &> {log};
        seqtk sample -2 -s {params.seed} {input.first_reads} \
            {params.number_to_sample} 2>> {log} | gzip > \
            {output.first_reads} 2>> {log};
        seqtk sample -2 -s {params.seed} {input.second_reads} \
            {params.number_to_sample} 2>> {log} | gzip > \
            {output.second_reads} 2>> {log}
        """
