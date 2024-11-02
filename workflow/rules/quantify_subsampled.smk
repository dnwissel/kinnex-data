configfile: "config/config.yaml"


singularity: f"docker://condaforge/mambaforge:{config['mambaforge_version']}"


rule quantify_subsampled_run_oarfish_subsampled:
    input:
        "results/downsample/run_minimap2_transcriptome_{type}/{subsample_number}/{read_number}/{sample}/{sample}.aligned.bam",
    output:
        "results/quantify_subsampled/run_oarfish_subsampled/{type}/{subsample_number}/{read_number}/{sample}/{sample}.quant",
    params:
        output_path="results/quantify_subsampled/run_oarfish_subsampled/{type}/{subsample_number}/{read_number}/{sample}/{sample}",
        filter_group=config["oarfish_filter_group"],
    log:
        "logs/quantify_subsampled/run_oarfish_subsampled/{type}/{subsample_number}/{read_number}/{sample}.out",
    threads: config["quantify_isoforms_threads"]
    conda:
        "../envs/standalone/oarfish.yaml"
    shell:
        """
        oarfish --version > {log};
        oarfish --threads {threads} \
                --filter-group {params.filter_group} \
                --model-coverage \
                --alignments {input} \
                --output {params.output_path} &>> {log}
        """


rule quantify_subsampled_run_salmon_illumina_subsampled:
    input:
        index="results/prepare/create_salmon_index/{type}",
        first_reads="results/subsample/convert_mapped_bams_to_fastq/{subsample_number}_{read_number}_{type}/{sample}-r1.fastq.gz",
        second_reads="results/subsample/convert_mapped_bams_to_fastq/{subsample_number}_{read_number}_{type}/{sample}-r2.fastq.gz",
    output:
        "results/quantify_subsampled/run_salmon_illumina_subsampled/{type}/{subsample_number}/{read_number}/{sample}/quant.sf",
    params:
        output_name="results/quantify_subsampled/run_salmon_illumina_subsampled/{type}/{subsample_number}/{read_number}/{sample}",
        fragment_length_mean=config["salmon_fragment_length_mean"],
        fragment_length_sd=config["salmon_fragment_length_sd"],
    log:
        "logs/quantify_subsampled/run_salmon_illumina_subsampled/{type}/{subsample_number}/{read_number}/{sample}.out",
    threads: config["quantify_isoforms_threads"]
    conda:
        "../envs/standalone/salmon.yaml"
    shell:
        """
        salmon --version > {log};
        salmon quant -i {input.index} -l A \
            -1 {input.first_reads} -2 {input.second_reads} \
            --validateMappings -o {params.output_name} -p {threads} \
            --seqBias --gcBias --fldMean {params.fragment_length_mean} \
            --fldSD {params.fragment_length_sd} &>> {log}
        """


rule quantify_subsampled_run_salmon_subsampled:
    input:
        reads="results/downsample/run_minimap2_transcriptome_{type}/{subsample_number}/{read_number}/{sample}/{sample}.aligned.bam",
        transcriptome="results/prepare/extract_transcriptomes/{type}_transcriptome.fa",
    output:
        "results/quantify_subsampled/run_salmon_subsampled/{type}/{subsample_number}/{read_number}/{sample}/quant.sf",
    params:
        output_path="results/quantify_subsampled/run_salmon_subsampled/{type}/{subsample_number}/{read_number}/{sample}",
        filter_group=config["oarfish_filter_group"],
    log:
        "logs/quantify_subsampled/run_salmon_subsampled/{type}/{subsample_number}/{read_number}/{sample}.out",
    threads: config["quantify_isoforms_threads"]
    conda:
        "../envs/standalone/salmon.yaml"
    shell:
        """
        salmon --version > {log};
        salmon quant --ont -t {input.transcriptome} -a {input.reads} \
                -o {params.output_path} -p {threads} -l A &>> {log}
        """


rule quantify_subsampled_run_bambu_subsampled:
    input:
        reads=expand(
            "results/downsample/run_minimap2_{{type}}/{{subsample_number}}/{{read_number}}/{sample}/{sample}.aligned.sorted.bam",
            sample=config["sample_names"],
        ),
        gencode_transcriptome="results/prepare/adjust_transcriptome_assembly_names/gencode.v45.primary_assembly.annotation.named.gtf",
        gencode_genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        sirv_transcriptome="results/prepare/adjust_sirv_names/sirv_set_four.gtf",
        sirv_genome="results/prepare/adjust_sirv_names/sirv_set_four.fa",
    output:
        "results/quantify_subsampled/run_bambu_subsampled/{type}/{subsample_number}/{read_number}/counts_transcript.txt",
    threads: config["quantify_isoforms_threads"]
    params:
        outdir="results/quantify_subsampled/run_bambu_subsampled/{type}/{subsample_number}/{read_number}",
        sirv="{type}",
    log:
        "logs/quantify_subsampled/run_bambu_subsampled/{type}/{subsample_number}/{read_number}/log.out",
    conda:
        "../envs/r/bambu.yaml"
    script:
        "../scripts/r/bambu_run_quantification.R"


rule quantify_subsampled_run_liqa_subsampled:
    input:
        reads=(
            "results/downsample/run_minimap2_sirv/{subsample_number}/{read_number}/{sample}/{sample}.aligned.sorted.bam",
        ),
        transcriptome="results/prepare/create_refgene_liqa_sirv/sirv.refgene",
    output:
        "results/quantify_subsampled/run_liqa_subsampled/sirv/{subsample_number}/{read_number}/{sample}.tsv",
    threads: config["quantify_isoforms_threads"]
    params:
        max_distance=config["liqa_max_distance"],
        f_weight=config["liqa_f_weight"],
    log:
        "logs/quantify_subsampled/run_liqa_subsampled/sirv/{subsample_number}/{read_number}/{sample}.out",
    conda:
        "../envs/standalone/liqa.yaml"
    shell:
        """
        liqa --version > {log};
        liqa -task quantify -refgene {input.transcriptome} -bam {input.reads} \
            -out {output} -max_distance {params.max_distance} \
            -f_weight {params.f_weight} &>> {log}
        """


rule quantify_subsampled_run_isoquant_subsampled:
    input:
        reads=expand(
            "results/downsample/run_minimap2_{{type}}/{{subsample_number}}/{{read_number}}/{sample}/{sample}.aligned.sorted.bam",
            sample=config["sample_names"],
        ),
        gencode_genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        sirv_transcriptome="results/prepare/make_db_files/sirv_set_four.db",
        gencode_transcriptome="results/prepare/make_db_files/gencode.v45.primary_assembly.annotation.named.db",
        sirv_genome="results/prepare/adjust_sirv_names/sirv_set_four.fa",
    output:
        directory(
            "results/quantify_subsampled/run_isoquant_subsampled/{type}/{subsample_number}/{read_number}"
        ),
    threads: config["quantify_isoforms_threads"]
    params:
        outdir="results/quantify_subsampled/run_isoquant_subsampled/{type}/{subsample_number}/{read_number}",
    log:
        "logs/quantify_subsampled/run_isoquant_subsampled/{type}/{subsample_number}/{read_number}/log.out",
    conda:
        "../envs/py/isoquant.yaml"
    shell:
        """
        isoquant.py --version > {log};
        if [ "{wildcards.type}" = "sirv" ]; then
            isoquant.py --reference {input.sirv_genome} \
                --genedb {input.sirv_transcriptome} \
                --bam {input.reads} \
                --polya_requirement never \
                --no_model_construction \
                --threads {threads} \
                --data_type pacbio_ccs -o {output} &>> {log}
        else
            isoquant.py --reference {input.gencode_genome} \
                --genedb {input.gencode_transcriptome} \
                --bam {input.reads} \
                --polya_requirement never \
                --no_model_construction \
                --threads {threads} \
                --data_type pacbio_ccs -o {output} &>> {log}
        fi
        """


rule quantify_subsampled_run_kallisto_long_subsampled:
    input:
        reads="results/subsample/convert_mapped_bams_to_fastq_long_reads/{subsample_number}_{read_number}_{type}/{sample}.fastq.gz",
        index=f"results/prepare/create_lr_kallisto_index/{{type}}_k-{config['lr_kallisto_index_k']}.idx",
        kallisto_binary="results/prepare/compile_lr_kallisto",
        gencode_transcriptome_gmap_headered="results/prepare/standardize_gtf_files/gencode_map_headered.txt",
        sirv_four_transcriptome_gmap_headered="results/prepare/standardize_gtf_files/sirv_set_four_map_headered.txt",
    output:
        directory(
            "results/quantify_subsampled/run_kallisto_long/{type}/{subsample_number}/{read_number}/{sample}"
        ),
    threads: config["quantify_isoforms_threads"]
    params:
        outdir="results/quantify_subsampled/run_kallisto_long/{type}/{subsample_number}/{read_number}/{sample}",
        bus_threshold=config["kallisto_bus_threshold"],
    log:
        "logs/quantify_subsampled/run_isoquant_subsampled/{type}/{subsample_number}/{read_number}/{sample}.out",
    conda:
        "../envs/standalone/kallisto.yaml"
    shell:
        """
        touch {log};
        if [ "{wildcards.type}" = "sirv" ]; then
            {input.kallisto_binary}/kallisto/build/src/kallisto  bus -t {threads} --long --threshold \
                {params.bus_threshold} -x bulk -i {input.index} \
                -o {params.outdir} {input.reads} &>> {log};
            bustools sort -t {threads} {params.outdir}/output.bus \
                -o {params.outdir}/sorted.bus &>> {log}; 
            bustools count {params.outdir}/sorted.bus \
            -t {params.outdir}/transcripts.txt \
            -e {params.outdir}/matrix.ec \
            -g {input.gencode_transcriptome_gmap_headered} \
            -o {params.outdir}/count --cm -m \
            &>> {log};
            {input.kallisto_binary}/kallisto/build/src/kallisto  quant-tcc -t {threads} \
                --long -P PacBio -f {params.outdir}/flens.txt \
                {params.outdir}/count.mtx -i {input.index} \
                -e {params.outdir}/count.ec.txt \
                -o {params.outdir} &>> {log}
        else
            {input.kallisto_binary}/kallisto/build/src/kallisto  bus -t {threads} --long --threshold \
                {params.bus_threshold} -x bulk -i {input.index} \
                -o {params.outdir} {input.reads} &>> {log};
            bustools sort -t {threads} {params.outdir}/output.bus \
                -o {params.outdir}/sorted.bus &>> {log}; 
            bustools count {params.outdir}/sorted.bus \
            -t {params.outdir}/transcripts.txt \
            -e {params.outdir}/matrix.ec \
            -g {input.gencode_transcriptome_gmap_headered} \
            -o {params.outdir}/count --cm -m \
            &>> {log};
            {input.kallisto_binary}/kallisto/build/src/kallisto  quant-tcc -t {threads} \
                --long -P PacBio -f {params.outdir}/flens.txt \
                {params.outdir}/count.mtx -i {input.index} \
                -e {params.outdir}/count.ec.txt \
                -o {params.outdir} &>> {log}
        fi
        """
