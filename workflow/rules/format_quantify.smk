configfile: "config/config.yaml"


singularity: f"docker://condaforge/mambaforge:{config['mambaforge_version']}"


SAMPLE_NAMES = [
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
]

SUBSAMPLED_SAMPLE_NAMES = [
    "day0-rep1",
    "day0-rep2",
    "day0-rep3",
    "day5-rep1",
    "day5-rep2",
    "day5-rep3",
]


rule format_quantify_isoquant:
    input:
        sample_paths=expand(
            "results/quantify/run_isoquant/{{type}}/{sample}/OUT/OUT.transcript_counts.tsv",
            sample=SAMPLE_NAMES,
        ),
        gencode_transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        sirv_transcriptome="results/prepare/standardize_gtf_files/sirv_set_four.gtf",
    output:
        output_path_transcript="results/format_quantify/isoquant/{type}/transcript_counts_formatted.tsv",
        output_path_gene="results/format_quantify/isoquant/{type}/gene_counts_formatted.tsv",
    params:
        transcript_id_col_ix=1,
        count_id_col_ix=2,
        sample_names=SAMPLE_NAMES,
    threads: 1
    log:
        "logs/format_quantify/isoquant/{type}.log",
    conda:
        "../envs/py/base.yaml"
    shell:
        """
        if [ "{wildcards.type}" = "sirv" ]; then
            python workflow/scripts/py/aggregate_counts.py \
                        --sample_paths {input.sample_paths} \
                        --sample_names {params.sample_names} \
                        --gtf_annotation_path {input.sirv_transcriptome} \
                        --output_path_transcript {output.output_path_transcript} \
                        --output_path_gene {output.output_path_gene} \
                        --transcript_id_col_ix {params.transcript_id_col_ix} \
                        --count_id_col_ix {params.count_id_col_ix}
        else
            python workflow/scripts/py/aggregate_counts.py \
                        --sample_paths {input.sample_paths} \
                        --sample_names {params.sample_names} \
                        --gtf_annotation_path {input.gencode_transcriptome} \
                        --output_path_transcript {output.output_path_transcript} \
                        --output_path_gene {output.output_path_gene} \
                        --transcript_id_col_ix {params.transcript_id_col_ix} \
                        --count_id_col_ix {params.count_id_col_ix}
        fi
        """


rule format_quantify_isoquant_novel:
    input:
        sample_paths=expand(
            "results/quantify/run_isoquant_novel/{{type}}/{sample}/OUT/OUT.transcript_counts.tsv",
            sample=SAMPLE_NAMES,
        ),
        gencode_novel_transcriptome="results/prepare/standardize_gtf_files_novel/gencode_novel.gtf",
    output:
        output_path_transcript="results/format_quantify/isoquant_novel/{type}/transcript_counts_formatted.tsv",
        output_path_gene="results/format_quantify/isoquant_novel/{type}/gene_counts_formatted.tsv",
    params:
        transcript_id_col_ix=1,
        count_id_col_ix=2,
        sample_names=SAMPLE_NAMES,
    threads: 1
    log:
        "logs/format_quantify/isoquant_novel/{type}.log",
    conda:
        "../envs/py/base.yaml"
    shell:
        """
        python workflow/scripts/py/aggregate_counts.py \
                    --sample_paths {input.sample_paths} \
                    --sample_names {params.sample_names} \
                    --gtf_annotation_path {input.gencode_novel_transcriptome} \
                    --output_path_transcript {output.output_path_transcript} \
                    --output_path_gene {output.output_path_gene} \
                    --transcript_id_col_ix {params.transcript_id_col_ix} \
                    --count_id_col_ix {params.count_id_col_ix}
        """


rule format_quantify_subsampled_isoquant:
    input:
        sample_paths=expand(
            "results/quantify_subsampled/run_isoquant/{{type}}/{{subsample_number}}/{{read_number}}/{sample}/OUT/OUT.transcript_counts.tsv",
            sample=SUBSAMPLED_SAMPLE_NAMES,
        ),
        gencode_transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        sirv_transcriptome="results/prepare/standardize_gtf_files/sirv_set_four.gtf",
    output:
        output_path_transcript="results/format_quantify_subsampled/isoquant/{type}/{subsample_number}/{read_number}/transcript_counts_formatted.tsv",
        output_path_gene="results/format_quantify_subsampled/isoquant/{type}/{subsample_number}/{read_number}/gene_counts_formatted.tsv",
    params:
        transcript_id_col_ix=1,
        count_id_col_ix=2,
        sample_names=SUBSAMPLED_SAMPLE_NAMES,
    threads: 1
    log:
        "logs/format_quantify_subsampled/isoquant/{type}/{subsample_number}/{read_number}.log",
    conda:
        "../envs/py/base.yaml"
    shell:
        """
        if [ "{wildcards.type}" = "sirv" ]; then
            python workflow/scripts/py/aggregate_counts.py \
                        --sample_paths {input.sample_paths} \
                        --sample_names {params.sample_names} \
                        --gtf_annotation_path {input.sirv_transcriptome} \
                        --output_path_transcript {output.output_path_transcript} \
                        --output_path_gene {output.output_path_gene} \
                        --transcript_id_col_ix {params.transcript_id_col_ix} \
                        --count_id_col_ix {params.count_id_col_ix}
        else
            python workflow/scripts/py/aggregate_counts.py \
                        --sample_paths {input.sample_paths} \
                        --sample_names {params.sample_names} \
                        --gtf_annotation_path {input.gencode_transcriptome} \
                        --output_path_transcript {output.output_path_transcript} \
                        --output_path_gene {output.output_path_gene} \
                        --transcript_id_col_ix {params.transcript_id_col_ix} \
                        --count_id_col_ix {params.count_id_col_ix}
        fi
        """


rule format_quantify_bambu:
    input:
        sample_paths=expand(
            "results/quantify/run_bambu/{{type}}/{sample}/counts_transcript.txt",
            sample=SAMPLE_NAMES,
        ),
        gencode_transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        sirv_transcriptome="results/prepare/standardize_gtf_files/sirv_set_four.gtf",
    output:
        output_path_transcript="results/format_quantify/bambu/{type}/transcript_counts_formatted.tsv",
        output_path_gene="results/format_quantify/bambu/{type}/gene_counts_formatted.tsv",
    params:
        transcript_id_col_ix=1,
        count_id_col_ix=3,
        sample_names=SAMPLE_NAMES,
    threads: 1
    log:
        "logs/format_quantify/bambu/{type}.log",
    conda:
        "../envs/py/base.yaml"
    shell:
        """
        if [ "{wildcards.type}" = "sirv" ]; then
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.sirv_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        else
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.gencode_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        fi
        """


rule format_quantify_bambu_novel:
    input:
        sample_paths=expand(
            "results/quantify/run_bambu_novel/{{type}}/{sample}/counts_transcript.txt",
            sample=SAMPLE_NAMES,
        ),
        gencode_novel_transcriptome="results/prepare/standardize_gtf_files_novel/gencode_novel.gtf",
    output:
        output_path_transcript="results/format_quantify/bambu_novel/{type}/transcript_counts_formatted.tsv",
        output_path_gene="results/format_quantify/bambu_novel/{type}/gene_counts_formatted.tsv",
    params:
        transcript_id_col_ix=1,
        count_id_col_ix=3,
        sample_names=SAMPLE_NAMES,
    threads: 1
    log:
        "logs/format_quantify/bambu_novel/{type}.log",
    conda:
        "../envs/py/base.yaml"
    shell:
        """
        python workflow/scripts/py/aggregate_counts.py \
                            --sample_paths {input.sample_paths} \
                            --sample_names {params.sample_names} \
                            --gtf_annotation_path {input.gencode_novel_transcriptome} \
                            --output_path_transcript {output.output_path_transcript} \
                            --output_path_gene {output.output_path_gene} \
                            --transcript_id_col_ix {params.transcript_id_col_ix} \
                            --count_id_col_ix {params.count_id_col_ix}
        """


rule format_quantify_subsampled_bambu:
    input:
        sample_paths=expand(
            "results/quantify_subsampled/run_bambu/{{type}}/{{subsample_number}}/{{read_number}}/{sample}/counts_transcript.txt",
            sample=SUBSAMPLED_SAMPLE_NAMES,
        ),
        gencode_transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        sirv_transcriptome="results/prepare/standardize_gtf_files/sirv_set_four.gtf",
    output:
        output_path_transcript="results/format_quantify_subsampled/bambu/{type}/{subsample_number}/{read_number}/transcript_counts_formatted.tsv",
        output_path_gene="results/format_quantify_subsampled/bambu/{type}/{subsample_number}/{read_number}/gene_counts_formatted.tsv",
    params:
        transcript_id_col_ix=1,
        count_id_col_ix=3,
        sample_names=SUBSAMPLED_SAMPLE_NAMES,
    threads: 1
    log:
        "logs/format_quantify_subsampled/bambu/{type}/{subsample_number}/{read_number}.log",
    conda:
        "../envs/py/base.yaml"
    shell:
        """
        if [ "{wildcards.type}" = "sirv" ]; then
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.sirv_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        else
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.gencode_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        fi
        """


rule format_quantify_oarfish:
    input:
        sample_paths=expand(
            "results/quantify/run_oarfish/{{type}}/{sample}/{sample}.quant",
            sample=SAMPLE_NAMES,
        ),
        gencode_transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        sirv_transcriptome="results/prepare/standardize_gtf_files/sirv_set_four.gtf",
    output:
        output_path_transcript="results/format_quantify/oarfish/{type}/transcript_counts_formatted.tsv",
        output_path_gene="results/format_quantify/oarfish/{type}/gene_counts_formatted.tsv",
    params:
        transcript_id_col_ix=1,
        count_id_col_ix=3,
        sample_names=SAMPLE_NAMES,
    threads: 1
    log:
        "logs/format_quantify/oarfish/{type}.log",
    conda:
        "../envs/py/base.yaml"
    shell:
        """
        if [ "{wildcards.type}" = "sirv" ]; then
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.sirv_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        else
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.gencode_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        fi
        """


rule format_quantify_oarfish_novel:
    input:
        sample_paths=expand(
            "results/quantify/run_oarfish_novel/{{type}}/{sample}/{sample}.quant",
            sample=SAMPLE_NAMES,
        ),
        gencode_novel_transcriptome="results/prepare/standardize_gtf_files_novel/gencode_novel.gtf",
    output:
        output_path_transcript="results/format_quantify/oarfish_novel/{type}/transcript_counts_formatted.tsv",
        output_path_gene="results/format_quantify/oarfish_novel/{type}/gene_counts_formatted.tsv",
    params:
        transcript_id_col_ix=1,
        count_id_col_ix=3,
        sample_names=SAMPLE_NAMES,
    threads: 1
    log:
        "logs/format_quantify/oarfish_novel/{type}.log",
    conda:
        "../envs/py/base.yaml"
    shell:
        """
        python workflow/scripts/py/aggregate_counts.py \
                            --sample_paths {input.sample_paths} \
                            --sample_names {params.sample_names} \
                            --gtf_annotation_path {input.gencode_novel_transcriptome} \
                            --output_path_transcript {output.output_path_transcript} \
                            --output_path_gene {output.output_path_gene} \
                            --transcript_id_col_ix {params.transcript_id_col_ix} \
                            --count_id_col_ix {params.count_id_col_ix}
        """


rule format_quantify_subsampled_oarfish:
    input:
        sample_paths=expand(
            "results/quantify_subsampled/run_oarfish/{{type}}/{{subsample_number}}/{{read_number}}/{sample}/{sample}.quant",
            sample=SUBSAMPLED_SAMPLE_NAMES,
        ),
        gencode_transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        sirv_transcriptome="results/prepare/standardize_gtf_files/sirv_set_four.gtf",
    output:
        output_path_transcript="results/format_quantify_subsampled/oarfish/{type}/{subsample_number}/{read_number}/transcript_counts_formatted.tsv",
        output_path_gene="results/format_quantify_subsampled/oarfish/{type}/{subsample_number}/{read_number}/gene_counts_formatted.tsv",
    params:
        transcript_id_col_ix=1,
        count_id_col_ix=3,
        sample_names=SUBSAMPLED_SAMPLE_NAMES,
    threads: 1
    log:
        "logs/format_quantify_subsampled/oarfish/{type}/{subsample_number}/{read_number}.log",
    conda:
        "../envs/py/base.yaml"
    shell:
        """
        if [ "{wildcards.type}" = "sirv" ]; then
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.sirv_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        else
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.gencode_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        fi
        """


rule format_quantify_salmon:
    input:
        sample_paths=expand(
            "results/quantify/run_salmon/{{type}}/{sample}/quant.sf",
            sample=SAMPLE_NAMES,
        ),
        gencode_transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        sirv_transcriptome="results/prepare/standardize_gtf_files/sirv_set_four.gtf",
    output:
        output_path_transcript="results/format_quantify/salmon/{type}/transcript_counts_formatted.tsv",
        output_path_gene="results/format_quantify/salmon/{type}/gene_counts_formatted.tsv",
    params:
        transcript_id_col_ix=1,
        count_id_col_ix=5,
        sample_names=SAMPLE_NAMES,
    threads: 1
    log:
        "logs/format_quantify/salmon/{type}.log",
    conda:
        "../envs/py/base.yaml"
    shell:
        """
        if [ "{wildcards.type}" = "sirv" ]; then
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.sirv_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        else
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.gencode_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        fi
        """


rule format_quantify_salmon_novel:
    input:
        sample_paths=expand(
            "results/quantify/run_salmon_novel/{{type}}/{sample}/quant.sf",
            sample=SAMPLE_NAMES,
        ),
        gencode_novel_transcriptome="results/prepare/standardize_gtf_files_novel/gencode_novel.gtf",
    output:
        output_path_transcript="results/format_quantify/salmon_novel/{type}/transcript_counts_formatted.tsv",
        output_path_gene="results/format_quantify/salmon_novel/{type}/gene_counts_formatted.tsv",
    params:
        transcript_id_col_ix=1,
        count_id_col_ix=5,
        sample_names=SAMPLE_NAMES,
    threads: 1
    log:
        "logs/format_quantify/salmon_novel/{type}.log",
    conda:
        "../envs/py/base.yaml"
    shell:
        """
        python workflow/scripts/py/aggregate_counts.py \
                            --sample_paths {input.sample_paths} \
                            --sample_names {params.sample_names} \
                            --gtf_annotation_path {input.gencode_novel_transcriptome} \
                            --output_path_transcript {output.output_path_transcript} \
                            --output_path_gene {output.output_path_gene} \
                            --transcript_id_col_ix {params.transcript_id_col_ix} \
                            --count_id_col_ix {params.count_id_col_ix}
        """


rule format_quantify_subsampled_salmon:
    input:
        sample_paths=expand(
            "results/quantify_subsampled/run_salmon/{{type}}/{{subsample_number}}/{{read_number}}/{sample}/quant.sf",
            sample=SUBSAMPLED_SAMPLE_NAMES,
        ),
        gencode_transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        sirv_transcriptome="results/prepare/standardize_gtf_files/sirv_set_four.gtf",
    output:
        output_path_transcript="results/format_quantify_subsampled/salmon/{type}/{subsample_number}/{read_number}/transcript_counts_formatted.tsv",
        output_path_gene="results/format_quantify_subsampled/salmon/{type}/{subsample_number}/{read_number}/gene_counts_formatted.tsv",
    params:
        transcript_id_col_ix=1,
        count_id_col_ix=5,
        sample_names=SUBSAMPLED_SAMPLE_NAMES,
    threads: 1
    log:
        "logs/format_quantify_subsampled/salmon/{type}/{subsample_number}/{read_number}.log",
    conda:
        "../envs/py/base.yaml"
    shell:
        """
        if [ "{wildcards.type}" = "sirv" ]; then
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.sirv_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        else
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.gencode_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        fi
        """


rule format_quantify_salmon_illumina:
    input:
        sample_paths=expand(
            "results/quantify/run_salmon_illumina/{{type}}/{sample}/quant.sf",
            sample=SAMPLE_NAMES,
        ),
        gencode_transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        sirv_transcriptome="results/prepare/standardize_gtf_files/sirv_set_four.gtf",
    output:
        output_path_transcript="results/format_quantify/salmon_illumina/{type}/transcript_counts_formatted.tsv",
        output_path_gene="results/format_quantify/salmon_illumina/{type}/gene_counts_formatted.tsv",
    params:
        transcript_id_col_ix=1,
        count_id_col_ix=5,
        sample_names=SAMPLE_NAMES,
    threads: 1
    log:
        "logs/format_quantify/salmon_illumina/{type}.log",
    conda:
        "../envs/py/base.yaml"
    shell:
        """
        if [ "{wildcards.type}" = "sirv" ]; then
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.sirv_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        else
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.gencode_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        fi
        """


rule format_quantify_salmon_illumina_novel:
    input:
        sample_paths=expand(
            "results/quantify/run_salmon_illumina_novel/{{type}}/{sample}/quant.sf",
            sample=SAMPLE_NAMES,
        ),
        gencode_novel_transcriptome="results/prepare/standardize_gtf_files_novel/gencode_novel.gtf",
    output:
        output_path_transcript="results/format_quantify/salmon_illumina_novel/{type}/transcript_counts_formatted.tsv",
        output_path_gene="results/format_quantify/salmon_illumina_novel/{type}/gene_counts_formatted.tsv",
    params:
        transcript_id_col_ix=1,
        count_id_col_ix=5,
        sample_names=SAMPLE_NAMES,
    threads: 1
    log:
        "logs/format_quantify/salmon_illumina_novel/{type}.log",
    conda:
        "../envs/py/base.yaml"
    shell:
        """
        python workflow/scripts/py/aggregate_counts.py \
                            --sample_paths {input.sample_paths} \
                            --sample_names {params.sample_names} \
                            --gtf_annotation_path {input.gencode_novel_transcriptome} \
                            --output_path_transcript {output.output_path_transcript} \
                            --output_path_gene {output.output_path_gene} \
                            --transcript_id_col_ix {params.transcript_id_col_ix} \
                            --count_id_col_ix {params.count_id_col_ix}
        """


rule format_quantify_subsampled_salmon_illumina:
    input:
        sample_paths=expand(
            "results/quantify_subsampled/run_salmon_illumina/{{type}}/{{subsample_number}}/{{read_number}}/{sample}/quant.sf",
            sample=SUBSAMPLED_SAMPLE_NAMES,
        ),
        gencode_transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        sirv_transcriptome="results/prepare/standardize_gtf_files/sirv_set_four.gtf",
    output:
        output_path_transcript="results/format_quantify_subsampled/salmon_illumina/{type}/{subsample_number}/{read_number}/transcript_counts_formatted.tsv",
        output_path_gene="results/format_quantify_subsampled/salmon_illumina/{type}/{subsample_number}/{read_number}/gene_counts_formatted.tsv",
    params:
        transcript_id_col_ix=1,
        count_id_col_ix=5,
        sample_names=SUBSAMPLED_SAMPLE_NAMES,
    threads: 1
    log:
        "logs/format_quantify_subsampled/salmon_illumina/{type}/{subsample_number}/{read_number}.log",
    conda:
        "../envs/py/base.yaml"
    shell:
        """
        if [ "{wildcards.type}" = "sirv" ]; then
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.sirv_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        else
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.gencode_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        fi
        """


rule format_quantify_kallisto:
    input:
        sample_paths=expand(
            "results/quantify/run_kallisto_long/{{type}}/{sample}/abundance_1.tsv",
            sample=SAMPLE_NAMES,
        ),
        gencode_transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        sirv_transcriptome="results/prepare/standardize_gtf_files/sirv_set_four.gtf",
    output:
        output_path_transcript="results/format_quantify/kallisto/{type}/transcript_counts_formatted.tsv",
        output_path_gene="results/format_quantify/kallisto/{type}/gene_counts_formatted.tsv",
    params:
        transcript_id_col_ix=1,
        count_id_col_ix=4,
        sample_names=SAMPLE_NAMES,
    threads: 1
    log:
        "logs/format_quantify/kallisto/{type}.log",
    conda:
        "../envs/py/base.yaml"
    shell:
        """
        if [ "{wildcards.type}" = "sirv" ]; then
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.sirv_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        else
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.gencode_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        fi
        """


rule format_quantify_kallisto_novel:
    input:
        sample_paths=expand(
            "results/quantify/run_kallisto_long_novel/{{type}}/{sample}/abundance_1.tsv",
            sample=SAMPLE_NAMES,
        ),
        gencode_novel_transcriptome="results/prepare/standardize_gtf_files_novel/gencode_novel.gtf",
    output:
        output_path_transcript="results/format_quantify/kallisto_novel/{type}/transcript_counts_formatted.tsv",
        output_path_gene="results/format_quantify/kallisto_novel/{type}/gene_counts_formatted.tsv",
    params:
        transcript_id_col_ix=1,
        count_id_col_ix=4,
        sample_names=SAMPLE_NAMES,
    threads: 1
    log:
        "logs/format_quantify/kallisto_novel/{type}.log",
    conda:
        "../envs/py/base.yaml"
    shell:
        """
        python workflow/scripts/py/aggregate_counts.py \
                            --sample_paths {input.sample_paths} \
                            --sample_names {params.sample_names} \
                            --gtf_annotation_path {input.gencode_novel_transcriptome} \
                            --output_path_transcript {output.output_path_transcript} \
                            --output_path_gene {output.output_path_gene} \
                            --transcript_id_col_ix {params.transcript_id_col_ix} \
                            --count_id_col_ix {params.count_id_col_ix}
        """


rule format_quantify_subsampled_kallisto:
    input:
        sample_paths=expand(
            "results/quantify_subsampled/run_kallisto_long/{{type}}/{{subsample_number}}/{{read_number}}/{sample}/abundance_1.tsv",
            sample=SUBSAMPLED_SAMPLE_NAMES,
        ),
        gencode_transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        sirv_transcriptome="results/prepare/standardize_gtf_files/sirv_set_four.gtf",
    output:
        output_path_transcript="results/format_quantify_subsampled/kallisto/{type}/{subsample_number}/{read_number}/transcript_counts_formatted.tsv",
        output_path_gene="results/format_quantify_subsampled/kallisto/{type}/{subsample_number}/{read_number}/gene_counts_formatted.tsv",
    params:
        transcript_id_col_ix=1,
        count_id_col_ix=4,
        sample_names=SUBSAMPLED_SAMPLE_NAMES,
    threads: 1
    log:
        "logs/format_quantify_subsampled/kallisto/{type}/{subsample_number}/{read_number}.log",
    conda:
        "../envs/py/base.yaml"
    shell:
        """
        if [ "{wildcards.type}" = "sirv" ]; then
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.sirv_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        else
            python workflow/scripts/py/aggregate_counts.py \
                                --sample_paths {input.sample_paths} \
                                --sample_names {params.sample_names} \
                                --gtf_annotation_path {input.gencode_transcriptome} \
                                --output_path_transcript {output.output_path_transcript} \
                                --output_path_gene {output.output_path_gene} \
                                --transcript_id_col_ix {params.transcript_id_col_ix} \
                                --count_id_col_ix {params.count_id_col_ix}
        fi
        """
