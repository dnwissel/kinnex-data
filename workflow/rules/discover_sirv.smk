configfile: "config/config.yaml"


rule discover_sirv_filter_annotation:
    input:
        "results/prepare/standardize_gtf_files/sirv_set_four.gtf",
    output:
        "results/discover_sirv/filter_annotation/sirv_set_four_filtered_short_sirvs.gtf",
    log:
        "logs/discover_sirv/filter_annotation/out.log",
    shell:
        """
        grep "SIRV[1-7]\\s" {input} > {output} 2> {log}
        """


rule discover_sirv_run_minimap2:
    input:
        reads="results/downsample/convert_mapped_bams_to_fastq_long_reads/{rep}_{num_reads}_sirv/{sample}.fastq.gz",
        sirv_genome="results/prepare/adjust_sirv_names/sirv_set_four.fa",
    output:
        "results/discover_sirv/run_minimap2/{rep}_{num_reads}/{sample}.aligned.sorted.bam",
    params:
        align_sort_bam_memory_gb=config["align_sort_bam_memory_gb"],
        align_sort_bam_threads=config["align_sort_bam_threads"],
    threads: config["align_map_bam_threads"]
    log:
        "logs/discover_sirv/run_minimap2/{rep}_{num_reads}/{sample}.log",
    conda:
        "../envs/standalone/minimap2.yaml"
    shell:
        """
        minimap2 -ax splice:hq \
            --splice-flank=no \
            -t {threads} -uf {input.sirv_genome} \
            {input.reads} 2>> {log} | \
            samtools sort -@ {params.align_sort_bam_threads} \
            -m{params.align_sort_bam_memory_gb}g -o {output} \
            - &>> {log}
        """


rule discover_sirv_run_bambu:
    input:
        reads=expand(
            "results/discover_sirv/run_minimap2/{{rep}}_{{num_reads}}/{sample}.aligned.sorted.bam",
            sample=["day0-rep1", "day0-rep2", "day0-rep3"],
        ),
        genome="results/prepare/adjust_sirv_names/sirv_set_four.fa",
    output:
        "results/discover_sirv/run_bambu_pb/{rep}_{num_reads}/transcriptome.gtf",
    log:
        "logs/discover_sirv/run_bambu/{rep}_{num_reads}/log.out",
    conda:
        "../envs/r/bambu.yaml"
    threads: config["discover_isoforms_threads"]
    script:
        "../scripts/r/discover/run_bambu_de_novo.R"


rule discover_sirv_index_star:
    input:
        genome="results/prepare/adjust_sirv_names/sirv_set_four.fa",
    output:
        directory("results/discover_sirv/index_star"),
    params:
        read_length=config["illumina_read_length"],
    threads: config["align_map_bam_threads"]
    log:
        "logs/discover_sirv/index_star/out.log",
    conda:
        "../envs/standalone/star.yaml"
    shell:
        """
        STAR --version > {log};
        STAR --runMode genomeGenerate --runThreadN {threads} \
            --genomeDir {output} --genomeFastaFiles {input.genome} \
            &>> {log}
        """


rule discover_sirv_run_illumina_star:
    input:
        reads_first="results/downsample/convert_mapped_bams_to_fastq/{rep}_{num_reads}_sirv/{sample}-r1.fastq.gz",
        reads_second="results/downsample/convert_mapped_bams_to_fastq/{rep}_{num_reads}_sirv/{sample}-r2.fastq.gz",
        genome="results/discover_sirv/index_star",
    output:
        "results/discover_sirv/run_illumina_star/{rep}_{num_reads}/{sample}_Aligned.sortedByCoord.out.bam",
    params:
        outfile_prefix="results/discover_sirv/run_illumina_star/{rep}_{num_reads}/{sample}_",
    threads: config["align_map_bam_threads"]
    log:
        "logs/discover_sirv/run_illumina_star/{rep}_{num_reads}/{sample}.log",
    conda:
        "../envs/standalone/star.yaml"
    shell:
        """
        STAR --version > {log};
        STAR --genomeDir {input.genome} --readFilesIn {input.reads_first} {input.reads_second} \
            --runThreadN {threads} --outFileNamePrefix {params.outfile_prefix} \
            --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif &>> {log};
        rm {input.reads_first} &>> {log};
        rm {input.reads_second} &>> {log}
        """


rule discover_sirv_index_illumina_star:
    input:
        "results/discover_sirv/run_illumina_star/{rep}_{num_reads}/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        "results/discover_sirv/run_illumina_star/{rep}_{num_reads}/{sample}_Aligned.sortedByCoord.out.bam.bai",
    log:
        "results/discover_sirv/index_illumina_star/{rep}_{num_reads}/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    shell:
        """
        samtools --version > {log};
        sleep 61s &>> {log};
        samtools index {input} >> {log}
        """


rule discover_sirv_run_stringtie2_illumina:
    input:
        reads="results/discover_sirv/run_illumina_star/{rep}_{num_reads}/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        "results/discover_sirv/run_stringtie2_illumina/{rep}_{num_reads}/{sample}.gtf",
    log:
        "logs/discover_sirv/run_stringtie2_illumina/{rep}_{num_reads}/{sample}.log",
    conda:
        "../envs/standalone/stringtie2.yaml"
    threads: config["discover_isoforms_threads"]
    shell:
        """
        conda list > {log};
        stringtie -p {threads} -o {output} {input.reads} &>> {log}
        """


rule discover_sirv_merge_stringtie2:
    input:
        discovered_transcriptomes=expand(
            "results/discover_sirv/run_stringtie2_illumina/{{rep}}_{{num_reads}}/{sample}.gtf",
            sample=["day0-rep1", "day0-rep2", "day0-rep3"],
        ),
    output:
        transcriptome="results/discover_sirv/merge_stringtie2_illumina/{rep}_{num_reads}/transcriptome.gtf",
        input_list="results/discover_sirv/merge_stringtie2_illumina/{rep}_{num_reads}/input_list.txt",
    log:
        "logs/discover_sirv/merge_stringtie2_illumina/{rep}_{num_reads}/out.log",
    conda:
        "../envs/standalone/stringtie2.yaml"
    shell:
        """
        conda list > {log};
        echo {input.discovered_transcriptomes} | tr " " "\n" > {output.input_list};
        stringtie --merge -i \
           -o {output.transcriptome} {output.input_list} \
           &>> {log}
        """


rule discover_sirv_evaluate_performance:
    input:
        transcriptome="results/discover_sirv/merge_stringtie2_illumina/{rep}_{num_reads}/transcriptome.gtf",
        ground_truth="results/discover_sirv/filter_annotation/sirv_set_four_filtered_short_sirvs.gtf",
    output:
        "results/discover_sirv/evaluate_performance/illumina_{rep}_{num_reads}",
    params:
        output_prefix="results/discover_sirv/evaluate_performance/illumina_{rep}_{num_reads}",
    log:
        "logs/discover_sirv/evaluate_performance/illumina_{rep}_{num_reads}.log",
    conda:
        "../envs/standalone/gffcompare.yaml"
    shell:
        """
        conda list > {log};
        gffcompare -r {input.ground_truth} -o {params.output_prefix} {input.transcriptome}
        """


rule discover_sirv_evaluate_performance_bambu:
    input:
        transcriptome="results/discover_sirv/run_bambu_pb/{rep}_{num_reads}/transcriptome.gtf",
        ground_truth="results/discover_sirv/filter_annotation/sirv_set_four_filtered_short_sirvs.gtf",
    output:
        "results/discover_sirv/evaluate_performance/pb_{rep}_{num_reads}",
    params:
        output_prefix="results/discover_sirv/evaluate_performance/pb_{rep}_{num_reads}",
    log:
        "logs/discover_sirv/evaluate_performance/pb_{rep}_{num_reads}.log",
    conda:
        "../envs/standalone/gffcompare.yaml"
    shell:
        """
        conda list > {log};
        gffcompare -r {input.ground_truth} -o {params.output_prefix} {input.transcriptome}
        """
