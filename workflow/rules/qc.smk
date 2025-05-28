configfile: "config/config.yaml"


singularity: f"docker://condaforge/mambaforge:{config['mambaforge_version']}"


rule qc_index_star_transcriptome:
    input:
        transcriptome="results/prepare/extract_transcriptomes/gencode_transcriptome.fa",
    output:
        directory("qc/index_star_transcriptome"),
    threads: config["align_map_bam_threads"]
    log:
        "logs/qc/index_star_transcriptome/out.log",
    conda:
        "../envs/standalone/star.yaml"
    shell:
        """
        STAR --version > {log};
        STAR --runMode genomeGenerate --runThreadN {threads} \
            --limitGenomeGenerateRAM 176852828426 \
            --genomeDir {output} --genomeFastaFiles {input.transcriptome} \
             &>> {log}
        """


rule qc_run_illumina_star_transcriptome:
    input:
        reads_first="results/align/convert_mapped_bams_to_fastq/{sample}/{data_type}/{sample}.reads_1.fastq.gz",
        reads_second="results/align/convert_mapped_bams_to_fastq/{sample}/{data_type}/{sample}.reads_2.fastq.gz",
        genome="results/prepare/index_star",
        transcriptome="results/prepare/concatenate_transcriptomes/transcriptome.gtf",
    output:
        "results/qc/run_illumina_star_transcriptome/{data_type}/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        "results/qc/run_illumina_star_transcriptome/{data_type}/{sample}/{sample}_Aligned.toTranscriptome.out.bam",
    params:
        outfile_prefix="results/qc/run_illumina_star_transcriptome/{data_type}/{sample}/{sample}_",
    threads: config["align_map_bam_threads"]
    log:
        "logs/qc/run_illumina_star_transcriptome/{data_type}/{sample}.out",
    conda:
        "../envs/standalone/star.yaml"
    shell:
        """
        STAR --version > {log};
        STAR --genomeDir {input.genome} --readFilesIn {input.reads_first} {input.reads_second} \
            --runThreadN {threads} --outFileNamePrefix {params.outfile_prefix} \
            --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c \
            --quantMode TranscriptomeSAM \
            --outSAMstrandField intronMotif &>> {log}
        """


rule qc_calculate_alignment_qc_illumina:
    input:
        "results/qc/sort_illumina_bam_name/{sample}/{type}/{sample}_Aligned.out.bam",
    output:
        "results/qc/calculate_alignment_qc_illumina/{sample}/{type}/alignment_qc.tsv",
    params:
        junction_regex="[2-9][0-9][0-9]*N",
    threads: 1
    log:
        "logs/qc/calculate_alignment_qc_illumina/{type}/{sample}.out",
    conda:
        "../envs/py/pysam.yaml"
    script:
        "../scripts/py/calculate_alignment_qc_ill.py"


rule qc_calculate_alignment_qc_kinnex:
    input:
        "results/align/run_minimap2_{type}/{sample}/{sample}.aligned.sorted.bam",
    output:
        "results/qc/calculate_alignment_qc_kinnex/{sample}/{type}/alignment_qc.tsv",
    params:
        junction_regex="[2-9][0-9][0-9]*N",
    threads: 12
    log:
        "logs/qc/calculate_alignment_qc_kinnex/{type}/{sample}.out",
    conda:
        "../envs/py/pysam.yaml"
    script:
        "../scripts/py/calculate_alignment_qc_kinnex.py"


rule qc_calculate_alignment_qc_illumina_transcriptome:
    input:
        "results/qc/sort_illumina_bam_name_transcriptome/{sample}/{type}/{sample}_Aligned.out.bam",
    output:
        "results/qc/calculate_alignment_qc_illumina_transcriptome/{sample}/{type}/alignment_qc.tsv",
    threads: 1
    log:
        "logs/qc/calculate_alignment_qc_illumina_transcriptome/{type}/{sample}.out",
    conda:
        "../envs/py/pysam.yaml"
    script:
        "../scripts/py/calculate_insert_size_ill.py"


rule qc_sort_illumina_bam_name:
    input:
        "results/qc/run_illumina_star_transcriptome/{data_type}/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        "results/qc/sort_illumina_bam_name/{sample}/{data_type}/{sample}_Aligned.out.bam",
    threads: 12
    log:
        "logs/qc/sort_illumina_bam_name/{data_type}/{sample}.out",
    conda:
        "../envs/standalone/samtools.yaml"
    shell:
        """
        samtools sort -n -@ 4 -m4g -o {output} {input}
        """


rule qc_sort_illumina_bam_name_transcriptome:
    input:
        "results/qc/run_illumina_star_transcriptome/{data_type}/{sample}/{sample}_Aligned.toTranscriptome.out.bam",
    output:
        "results/qc/sort_illumina_bam_name_transcriptome/{sample}/{data_type}/{sample}_Aligned.out.bam",
    threads: 12
    log:
        "logs/qc/sort_illumina_bam_name_transcriptome/{data_type}/{sample}.out",
    conda:
        "../envs/standalone/samtools.yaml"
    shell:
        """
        samtools sort -n -@ 4 -m4g -o {output} {input}
        """


rule qc_calculate_alignment_qc_kinnex_transcriptome:
    input:
        "results/align/run_minimap2_transcriptome_{type}/{sample}/{sample}.aligned.bam",
    output:
        "results/qc/calculate_alignment_qc_kinnex_transcriptome/{sample}/{type}/alignment_qc.tsv",
    threads: 12
    log:
        "logs/qc/calculate_alignment_qc_kinnex_transcriptome/{type}/{sample}.out",
    conda:
        "../envs/py/pysam.yaml"
    script:
        "../scripts/py/calculate_insert_size_kinnex.py"


rule qc_prepare_sampled_read_quality_frame_illumina:
    input:
        input_path_transcriptome="results/qc/calculate_alignment_qc_illumina_transcriptome/{sample}/{type}/alignment_qc.tsv",
        first_lengths="results/qc/calculate_read_lengths_illumina/{sample}/{type}/lengths_1.txt",
        second_lengths="results/qc/calculate_read_lengths_illumina/{sample}/{type}/lengths_2.txt",
        first_quality="results/qc/calculate_base_q_illumina/{sample}/{type}/quality_1.txt",
        second_quality="results/qc/calculate_base_q_illumina/{sample}/{type}/quality_2.txt",
        alignment="results/qc/calculate_alignment_qc_illumina/{sample}/{type}/alignment_qc.tsv",
    output:
        "results/qc/prepare_sampled_read_quality_frame_illumina/{sample}/{type}/quality_plot_frame.tsv",
    params:
        seed=config["seed"],
        sample_number=config["qc_sample_number"],
    threads: 12
    log:
        "logs/qc/prepare_sampled_read_quality_frame_illumina/{type}/{sample}.out",
    conda:
        "../envs/r/qual_qc.yaml"
    script:
        "../scripts/r/qc_prepare_sampled_read_quality_frame_illumina.R"


rule qc_prepare_sampled_read_quality_frame_kinnex:
    input:
        lengths="results/qc/calculate_read_lengths_pb/{sample}/{type}/lengths.txt",
        quality="results/qc/calculate_base_q_pb/{sample}/{type}/quality.txt",
        alignment="results/qc/calculate_alignment_qc_kinnex/{sample}/{type}/alignment_qc.tsv",
    output:
        "results/qc/prepare_sampled_read_quality_frame_kinnex/{sample}/{type}/quality_plot_frame.tsv",
    params:
        seed=config["seed"],
        sample_number=config["qc_sample_number"],
    threads: 12
    log:
        "logs/qc/prepare_sampled_read_quality_frame_kinnex/{type}/{sample}.out",
    conda:
        "../envs/r/qual_qc.yaml"
    script:
        "../scripts/r/qc_prepare_sampled_read_quality_frame_kinnex.R"


rule qc_calculate_read_lengths_illumina:
    input:
        first_reads="results/align/convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.reads_1.fastq.gz",
        second_reads="results/align/convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.reads_2.fastq.gz",
    output:
        first_lengths="results/qc/calculate_read_lengths_illumina/{sample}/{type}/lengths_1.txt",
        second_lengths="results/qc/calculate_read_lengths_illumina/{sample}/{type}/lengths_2.txt",
    threads: 12
    log:
        "logs/qc/calculate_read_lengths_illumina/{type}/{sample}.out",
    conda:
        "../envs/standalone/qc.yaml"
    shell:
        """
        conda list > {log};
        bioawk -c fastx '{{ print $name, length($seq) }}' {input.first_reads} > \
            {output.first_lengths} 2>> {log};
        bioawk -c fastx '{{ print $name, length($seq) }}' {input.second_reads} > \
            {output.second_lengths} 2>> {log}
        """


rule qc_calculate_read_lengths_pb:
    input:
        "results/prepare/convert_mapped_bams_to_fastq/{sample}/{type}.reads.fastq.gz",
    output:
        "results/qc/calculate_read_lengths_pb/{sample}/{type}/lengths.txt",
    threads: 12
    log:
        "logs/qc/calculate_read_lengths_pb/{type}/{sample}.out",
    conda:
        "../envs/standalone/qc.yaml"
    shell:
        """
        conda list > {log};
        bioawk -c fastx '{{ print $name, length($seq) }}' {input} > \
            {output} 2>> {log}
        """


rule qc_calculate_base_q_illumina:
    input:
        first_reads="results/align/convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.reads_1.fastq.gz",
        second_reads="results/align/convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.reads_2.fastq.gz",
    output:
        first_quality="results/qc/calculate_base_q_illumina/{sample}/{type}/quality_1.txt",
        second_quality="results/qc/calculate_base_q_illumina/{sample}/{type}/quality_2.txt",
    threads: 12
    log:
        "logs/qc/calculate_base_q_illumina/{type}/{sample}.out",
    conda:
        "../envs/standalone/qc.yaml"
    shell:
        """
        conda list > {log};
        bioawk -c fastx '{{ print $name, meanqual($qual) }}' \
            {input.first_reads} > {output.first_quality} 2>> {log};
        bioawk -c fastx '{{ print $name, meanqual($qual) }}' \
            {input.second_reads} > {output.second_quality} 2>> {log}
        """


rule qc_calculate_base_q_pb:
    input:
        "results/prepare/convert_mapped_bams_to_fastq/{sample}/{type}.reads.fastq.gz",
    output:
        "results/qc/calculate_base_q_pb/{sample}/{type}/quality.txt",
    threads: 12
    log:
        "logs/qc/calculate_base_q_pb/{type}/{sample}.out",
    conda:
        "../envs/standalone/qc.yaml"
    shell:
        """
        conda list > {log};
        bioawk -c fastx '{{ print $name, meanqual($qual) }}' \
            {input} > {output} 2>> {log}
        """


rule qc_subsample_transcripts_three_prime_bias:
    input:
        transcriptome="results/prepare/adjust_transcriptome_assembly_names/gencode.v45.primary_assembly.annotation.named.gtf",
    output:
        subsampled_transcriptome="results/qc/subsample_transcripts_three_prime_bias/gencode_subsampled.gtf",
    params:
        n_subsample=config["three_prime_bias_subsample"],
        seed=config["seed"],
    log:
        "logs/qc/subsample_transcripts_three_prime_bias/log.out",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/subsample_transcripts_three_prime_bias.R"


rule qc_convert_gtf_to_bed:
    input:
        gencode_transcriptome="results/qc/subsample_transcripts_three_prime_bias/gencode_subsampled.gtf",
        sirv_transcriptome="results/prepare/adjust_sirv_names/sirv_set_four.gtf",
    output:
        gencode_bed="results/qc/convert_gtf_to_bed/gencode.bed",
        sirv_bed="results/qc/convert_gtf_to_bed/sirv.bed",
    log:
        "logs/qc/convert_gtf_to_bed/log.out",
    conda:
        "../envs/standalone/ucsc_genepred.yaml"
    shell:
        """
        conda list &> {log};
        cat {input.gencode_transcriptome} 2>> {log} |\
            gtfToGenePred /dev/stdin /dev/stdout 2>> {log} |\
            genePredToBed /dev/stdin /dev/stdout > \
            {output.gencode_bed} 2>> {log};
        cat {input.sirv_transcriptome} 2>> {log} |\
            gtfToGenePred /dev/stdin /dev/stdout 2>> {log} |\
            genePredToBed /dev/stdin /dev/stdout > \
            {output.sirv_bed} 2>> {log}
        """


rule qc_calculate_three_prime_bias_pb:
    input:
        gencode_bed="results/qc/convert_gtf_to_bed/gencode.bed",
        sirv_bed="results/qc/convert_gtf_to_bed/sirv.bed",
        bam="results/align/run_minimap2_{type}/{sample}/{sample}.aligned.sorted.bam",
    output:
        "results/qc/calculate_three_prime_bias_pb/{type}/{sample}/three_prime_bias.geneBodyCoverage.txt",
    params:
        output_path="results/qc/calculate_three_prime_bias_pb/{type}/{sample}/three_prime_bias",
    threads: 12
    log:
        "logs/qc/calculate_three_prime_bias_pb/{type}/{sample}.out",
    conda:
        "../envs/standalone/rseqc.yaml"
    shell:
        """
        conda list > {log};
        if [ "{wildcards.type}" = "sirv" ]; then
            geneBody_coverage.py \
                -r {input.sirv_bed} \
                -o {params.output_path} \
                -i {input.bam} &>> {log}
        else
            geneBody_coverage.py \
                -r {input.gencode_bed} \
                -o {params.output_path} \
                -i {input.bam} &>> {log}
        fi
        """


rule qc_calculate_three_prime_bias_illumina:
    input:
        gencode_bed="results/qc/convert_gtf_to_bed/gencode.bed",
        sirv_bed="results/qc/convert_gtf_to_bed/sirv.bed",
        bam="results/align/separate_sirv_and_gencode_illumina/{sample}/{sample}.aligned.{type}.sorted.bam",
    output:
        "results/qc/calculate_three_prime_bias_illumina/{type}/{sample}/three_prime_bias.geneBodyCoverage.txt",
    params:
        output_path="results/qc/calculate_three_prime_bias_illumina/{type}/{sample}/three_prime_bias",
    log:
        "logs/qc/calculate_three_prime_bias_illumina/{type}/{sample}.out",
    conda:
        "../envs/standalone/rseqc.yaml"
    shell:
        """
        conda list > {log};
        if [ "{wildcards.type}" = "sirv" ]; then
            geneBody_coverage.py \
                -r {input.sirv_bed} \
                -o {params.output_path} \
                -i {input.bam} &>> {log}
        else
            geneBody_coverage.py \
                -r {input.gencode_bed} \
                -o {params.output_path} \
                -i {input.bam} &>> {log}
        fi
        """
