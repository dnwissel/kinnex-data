configfile: "config/config.yaml"


singularity: f"docker://condaforge/mambaforge:{config['mambaforge_version']}"


rule qc_calculate_read_lengths_illumina:
    input:
        first_reads="results/align/convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.reads_1.fastq.gz",
        second_reads="results/align/convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.reads_2.fastq.gz",
    output:
        first_lengths="results/qc/calculate_read_lengths_illumina/{sample}/{type}/lengths_1.txt",
        second_lengths="results/qc/calculate_read_lengths_illumina/{sample}/{type}/lengths_2.txt",
    log:
        "logs/qc/calculate_read_lengths_illumina/{type}/{sample}.out",
    conda:
        "../envs/standalone/qc.yaml"
    shell:
        """
        conda list > {log};
        zcat {input.first_reads} | \
            awk '{{if(NR%4==2) print length($1)}}' > \
            {output.first_lengths} 2>> {log};
        zcat {input.second_reads} | \
            awk '{{if(NR%4==2) print length($1)}}' > \
            {output.second_lengths} 2>> {log}
        """


rule qc_calculate_read_lengths_pb:
    input:
        "results/prepare/convert_mapped_bams_to_fastq/{sample}/{type}.reads.fastq.gz",
    output:
        "results/qc/calculate_read_lengths_pb/{sample}/{type}/lengths.txt",
    log:
        "logs/qc/calculate_read_lengths_pb/{type}/{sample}.out",
    conda:
        "../envs/standalone/qc.yaml"
    shell:
        """
        conda list > {log};
        zcat {input} | \
            awk '{{if(NR%4==2) print length($1)}}' > \
            {output} 2>> {log}
        """


rule qc_calculate_base_q_illumina:
    input:
        first_reads="results/align/convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.reads_1.fastq.gz",
        second_reads="results/align/convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.reads_2.fastq.gz",
    output:
        first_quality="results/qc/calculate_base_q_illumina/{sample}/{type}/quality_1.fasta",
        second_quality="results/qc/calculate_base_q_illumina/{sample}/{type}/quality_2.fasta",
    log:
        "logs/qc/calculate_base_q_illumina/{type}/{sample}.out",
    conda:
        "../envs/standalone/qc.yaml"
    shell:
        """
        conda list > {log};
        bioawk -c fastx '{{print ">"$name; print meanqual($qual)}}' \
            {input.first_reads} > {output.first_quality} 2>> {log};
        bioawk -c fastx '{{print ">"$name; print meanqual($qual)}}' \
            {input.second_reads} > {output.second_quality} 2>> {log}
        """


rule qc_calculate_base_q_pb:
    input:
        "results/prepare/convert_mapped_bams_to_fastq/{sample}/{type}.reads.fastq.gz",
    output:
        "results/qc/calculate_base_q_pb/{sample}/{type}/quality.txt",
    log:
        "logs/qc/calculate_base_q_pb/{type}/{sample}.out",
    conda:
        "../envs/standalone/qc.yaml"
    shell:
        """
        conda list > {log};
        bioawk -c fastx '{{print ">"$name; print meanqual($qual)}}' \
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
