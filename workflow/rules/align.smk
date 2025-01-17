configfile: "config/config.yaml"


rule align_run_minimap2:
    input:
        reads="results/prepare/convert_mapped_bams_to_fastq/{sample}/{type}.reads.fastq.gz",
        gencode_genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        sirv_genome="results/prepare/adjust_sirv_names/sirv_set_four.fa",
        transcriptome="results/prepare/convert_gtfs_to_beds/{type}.bed",
    output:
        "results/align/run_minimap2_{type}/{sample}/{sample}.aligned.sorted.bam",
    params:
        align_sort_bam_memory_gb=config["align_sort_bam_memory_gb"],
        align_sort_bam_threads=config["align_sort_bam_threads"],
    threads: config["align_map_bam_threads"]
    log:
        "logs/align/run_minimap2/{type}/{sample}.log",
    conda:
        "../envs/standalone/minimap2.yaml"
    shell:
        """
        minimap2 --version > {log};
        if [ "{wildcards.type}" = "sirv" ]; then
            minimap2 -ax splice:hq --junc-bed {input.transcriptome} \
                --splice-flank=no \
                -t {threads} -uf {input.sirv_genome} \
                {input.reads} 2>> {log} | \
                samtools sort -@ {params.align_sort_bam_threads} \
                -m{params.align_sort_bam_memory_gb}g -o {output} \
                - &>> {log};
            sleep 61s &>> {log};
            samtools index {output} &>> {log}
        else
            minimap2 -ax splice:hq --junc-bed {input.transcriptome} \
                -t {threads} -uf {input.gencode_genome} \
                {input.reads} 2>> {log} | \
                samtools sort -@ {params.align_sort_bam_threads} \
                -m{params.align_sort_bam_memory_gb}g -o {output} \
                - &>> {log};
            sleep 61s &>> {log};
            samtools index {output} &>> {log}
        fi
        """


# rule align_run_minimap2_transcriptome:
#     input:
#         reads="results/prepare/convert_mapped_bams_to_fastq/{sample}/{type}.reads.fastq.gz",
#         transcriptome="results/prepare/extract_transcriptomes/{type}_transcriptome.fa",
#     output:
#         "results/align/run_minimap2_transcriptome_{type}/{sample}/{sample}.aligned.bam",
#     params:
#         n_secondary_alignments=config["n_secondary_alignments"],
#         memory=config["align_sort_bam_memory_gb"],
#         align_sort_bam_threads=config["align_sort_bam_threads"],
#     threads: config["align_map_bam_threads"]
#     log:
#         "logs/align/run_minimap2_transcriptome_gencode/{type}/{sample}.log",
#     benchmark:
#         repeat(
#             "benchmarks/align/run_minimap2_transcriptome_gencode/{type}/{sample}.txt",
#             config["timing_repetitions"],
#         )
#     conda:
#         "../envs/standalone/minimap2.yaml"
#     shell:
#         """
#         minimap2 --version > {log};
#         minimap2 --eqx -N {params.n_secondary_alignments} -ax map-hifi \
#             -t {threads} {input.transcriptome} \
#             {input.reads} 2>> {log} | \
#             samtools sort -n -@ {params.align_sort_bam_threads} \
#             -m{params.memory}g -o {output} - &>> {log}
#         """


rule align_run_illumina_star:
    input:
        reads_first="results/prepare/run_trim_galore/{sample}-r1_val_1.fq.gz",
        reads_second="results/prepare/run_trim_galore/{sample}-r2_val_2.fq.gz",
        genome="results/prepare/index_star",
    output:
        "results/align/run_illumina_star/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    params:
        outfile_prefix="results/align/run_illumina_star/{sample}/{sample}_",
    threads: config["align_map_bam_threads"]
    log:
        "logs/align/run_illumina_star/{sample}.out",
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


rule align_index_illumina_star:
    input:
        "results/align/run_illumina_star/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        "results/align/run_illumina_star/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/align/index_illumina_star/{sample}.out",
    conda:
        "../envs/standalone/samtools.yaml"
    shell:
        """
        samtools --version > {log};
        sleep 61s &>> {log};
        samtools index {input} >> {log}
        """


rule align_separate_sirv_and_gencode_illumina:
    input:
        bam="results/align/run_illumina_star/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        index="results/align/run_illumina_star/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai",
    output:
        gencode_reads="results/align/separate_sirv_and_gencode_illumina/{sample}/{sample}.aligned.gencode.sorted.bam",
        gencode_ix="results/align/separate_sirv_and_gencode_illumina/{sample}/{sample}.aligned.gencode.sorted.bam.bai",
        sirv_reads="results/align/separate_sirv_and_gencode_illumina/{sample}/{sample}.aligned.sirv.sorted.bam",
        sirv_ix="results/align/separate_sirv_and_gencode_illumina/{sample}/{sample}.aligned.sirv.sorted.bam.bai",
    params:
        lines_to_cut=1,
        gencode_chromosomes="chr",
    threads: config["stall_io_threads"]
    log:
        "logs/align/separate_sirv_and_gencode_illumina/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    shell:
        """
        samtools --version > {log};
        samtools idxstats {input.bam} 2>> {log} | \
            cut -f {params.lines_to_cut} 2>> {log} | \
            grep {params.gencode_chromosomes}  2>> {log} | \
            xargs samtools view -F 0x904 -b {input.bam} | samtools view -b -f 0x2 - > {output.gencode_reads} 2>> {log};
        samtools idxstats {input.bam} 2>> {log} | \
            cut -f {params.lines_to_cut} 2>> {log} | \
            grep -v {params.gencode_chromosomes} 2>> {log} | \
            head -n -{params.lines_to_cut} 2>> {log} | \
            xargs samtools view -F 0x904 -b {input.bam} | samtools view -b -f 0x2 - > {output.sirv_reads} 2>> {log};
        sleep 61s &>> {log};
        samtools index {output.gencode_reads} &>> {log};
        samtools index {output.sirv_reads} &>> {log};
        rm {input.bam} &>> {log};
        rm {input.index} &>> {log}
        """


rule align_convert_illumina_mapped_bams_to_fastq:
    input:
        "results/align/separate_sirv_and_gencode_illumina/{sample}/{sample}.aligned.{type}.sorted.bam",
    output:
        singletons="results/align/convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.singletons.fastq.gz",
        first_reads="results/align/convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.reads_1.fastq.gz",
        second_reads="results/align/convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.reads_2.fastq.gz",
    params:
        memory=config["align_sort_bam_memory_gb"],
        compression_level=config["compression_level"],
    threads: config["stall_io_threads"]
    log:
        "logs/prepare/filter_genome_mappings/{type}/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    shell:
        """
        samtools --version > {log};
        samtools sort -n -m{params.memory}g -@ {threads} {input} 2>> {log} | \
            samtools fastq -@ {threads} -1 {output.first_reads} \
            -2 {output.second_reads} -s {output.singletons} -c {params.compression_level} \
            - &>> {log}
        """
