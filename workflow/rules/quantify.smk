configfile: "config/config.yaml"


singularity: f"docker://condaforge/mambaforge:{config['mambaforge_version']}"


rule quantify_run_oarfish:
    input:
        "results/align/run_minimap2_transcriptome_{type}/{sample}/{sample}.aligned.bam",
    output:
        "results/quantify/run_oarfish/{type}/{sample}/{sample}.quant",
    params:
        output_path="results/quantify/run_oarfish/{type}/{sample}/{sample}",
        filter_group=config["oarfish_filter_group"],
    threads: config["quantify_isoforms_threads"]
    log:
        "logs/quantify/run_oarfish/{type}/{sample}.out",
    benchmark:
        repeat(
            "benchmarks/quantify/run_oarfish/{type}/{sample}.txt",
            config["timing_repetitions"],
        )
    conda:
        "../envs/standalone/oarfish.yaml"
    shell:
        """
        conda list > {log};
        oarfish --threads {threads} \
                --filter-group {params.filter_group} \
                --model-coverage \
                --alignments {input} \
                --output {params.output_path} &>> {log}
        """


rule quantify_run_salmon:
    input:
        reads="results/align/run_minimap2_transcriptome_{type}/{sample}/{sample}.aligned.bam",
        transcriptome="results/prepare/extract_transcriptomes/{type}_transcriptome.fa",
    output:
        "results/quantify/run_salmon/{type}/{sample}/quant.sf",
    params:
        output_path="results/quantify/run_salmon/{type}/{sample}",
    threads: config["quantify_isoforms_threads"]
    log:
        "logs/quantify/run_salmon/{type}/{sample}.out",
    benchmark:
        repeat(
            "benchmarks/quantify/run_salmon/{type}/{sample}.txt",
            config["timing_repetitions"],
        )
    conda:
        "../envs/standalone/salmon.yaml"
    shell:
        """
        salmon --version > {log};
        salmon quant --ont -t {input.transcriptome} -a {input.reads} \
                -o {params.output_path} -p {threads} -l A &>> {log}
        """


rule quantify_run_bambu:
    input:
        reads="results/align/run_minimap2_{type}/{sample}/{sample}.aligned.sorted.bam",
        gencode_transcriptome="results/prepare/adjust_transcriptome_assembly_names/gencode.v45.primary_assembly.annotation.named.gtf",
        gencode_genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        sirv_transcriptome="results/prepare/adjust_sirv_names/sirv_set_four.gtf",
        sirv_genome="results/prepare/adjust_sirv_names/sirv_set_four.fa",
    output:
        "results/quantify/run_bambu/{type}/{sample}/counts_transcript.txt",
    threads: config["quantify_isoforms_threads"]
    params:
        outdir="results/quantify/run_bambu/{type}/{sample}",
        sirv="{type}",
    log:
        "logs/quantify/run_bambu/{type}/{sample}.out",
    benchmark:
        repeat(
            "benchmarks/quantify/run_bambu/{type}/{sample}.txt",
            config["timing_repetitions"],
        )
    conda:
        "../envs/r/bambu.yaml"
    script:
        "../scripts/r/bambu_run_quantification.R"


rule quantify_run_isoquant:
    input:
        reads="results/align/run_minimap2_{type}/{sample}/{sample}.aligned.sorted.bam",
        gencode_genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        sirv_transcriptome="results/prepare/make_db_files/sirv_set_four.db",
        gencode_transcriptome="results/prepare/make_db_files/gencode.v45.primary_assembly.annotation.named.db",
        sirv_genome="results/prepare/adjust_sirv_names/sirv_set_four.fa",
    output:
        directory("results/quantify/run_isoquant/{type}/{sample}"),
    threads: config["quantify_isoforms_threads"]
    params:
        outdir="results/quantify/run_isoquant/{type}/{sample}",
    log:
        "logs/quantify/run_isoquant/{type}/{sample}.out",
    benchmark:
        repeat(
            "benchmarks/quantify/run_isoquant/{type}/{sample}.txt",
            config["timing_repetitions"],
        )
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
                --complete_genedb \
                --threads {threads} \
                --data_type pacbio_ccs -o {output} &>> {log}
        else
            isoquant.py --reference {input.gencode_genome} \
                --genedb {input.gencode_transcriptome} \
                --bam {input.reads} \
                --polya_requirement never \
                --no_model_construction \
                --complete_genedb \
                --threads {threads} \
                --data_type pacbio_ccs -o {output} &>> {log}
        fi
        """


rule quantify_run_kallisto_long:
    input:
        reads="results/prepare/convert_mapped_bams_to_fastq/{sample}/{type}.reads.fastq.gz",
        index=f"results/prepare/create_lr_kallisto_index/{{type}}_k-{config['lr_kallisto_index_k']}.idx",
        kallisto_binary="results/prepare/compile_lr_kallisto",
        gencode_transcriptome_gmap_headered="results/prepare/standardize_gtf_files/gencode_map_headered.txt",
        sirv_four_transcriptome_gmap_headered="results/prepare/standardize_gtf_files/sirv_set_four_map_headered.txt",
    output:
        directory("results/quantify/run_kallisto_long/{type}/{sample}"),
    threads: config["quantify_isoforms_threads"]
    params:
        outdir="results/quantify/run_kallisto_long/{type}/{sample}",
        bus_threshold=config["kallisto_bus_threshold"],
    log:
        "logs/quantify/run_kallisto_long/{type}/{sample}.out",
    benchmark:
        repeat(
            "benchmarks/quantify/run_kallisto_long/{type}/{sample}.txt",
            config["timing_repetitions"],
        )
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


rule quantify_run_salmon_illumina:
    input:
        index="results/prepare/create_salmon_index/{type}",
        first_reads="results/align/convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.reads_1.fastq.gz",
        second_reads="results/align/convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.reads_2.fastq.gz",
    output:
        "results/quantify/run_salmon_illumina/{type}/{sample}/quant.sf",
    params:
        output_name="results/quantify/run_salmon_illumina/{type}/{sample}",
        fragment_length_mean=config["salmon_fragment_length_mean"],
        fragment_length_sd=config["salmon_fragment_length_sd"],
    log:
        "logs/quantify/run_salmon_illumina/{type}/{sample}.out",
    threads: config["quantify_isoforms_threads"]
    conda:
        "../envs/standalone/salmon.yaml"
    shell:
        """
        conda list > {log};
        salmon quant -i {input.index} -l A \
            -1 {input.first_reads} -2 {input.second_reads} \
            --validateMappings -o {params.output_name} -p {threads} \
            --seqBias --gcBias --fldMean {params.fragment_length_mean} \
            --fldSD {params.fragment_length_sd} &>> {log}
        """


rule quantify_run_salmon_illumina_novel:
    input:
        index="results/prepare/create_salmon_index_novel/{type}",
        first_reads="results/align/convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.reads_1.fastq.gz",
        second_reads="results/align/convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.reads_2.fastq.gz",
    output:
        "results/quantify/run_salmon_illumina_novel/{type}/{sample}/quant.sf",
    params:
        output_name="results/quantify/run_salmon_illumina_novel/{type}/{sample}",
        fragment_length_mean=config["salmon_fragment_length_mean"],
        fragment_length_sd=config["salmon_fragment_length_sd"],
    log:
        "logs/quantify/run_salmon_illumina_novel/{type}/{sample}.out",
    threads: config["quantify_isoforms_threads"]
    conda:
        "../envs/standalone/salmon.yaml"
    shell:
        """
        conda list > {log};
        salmon quant -i {input.index} -l A \
            -1 {input.first_reads} -2 {input.second_reads} \
            --validateMappings -o {params.output_name} -p {threads} \
            --seqBias --gcBias --fldMean {params.fragment_length_mean} \
            --fldSD {params.fragment_length_sd} &>> {log}
        """
