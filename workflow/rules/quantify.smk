configfile: "config/config.yaml"


singularity: f"docker://condaforge/mambaforge:{config['mambaforge_version']}"


wildcard_constraints:
    data_type="(sirv|gencode)",
    sample="day[0-5]-rep[1-3]",


rule quantify_run_oarfish:
    input:
        "results/align/run_minimap2_transcriptome_{data_type}/{sample}/{sample}.aligned.bam",
    output:
        "results/quantify/run_oarfish/{data_type}/{sample}/{sample}.quant",
    params:
        output_path="results/quantify/run_oarfish/{data_type}/{sample}/{sample}",
        filter_group=config["oarfish_filter_group"],
    threads: config["quantify_isoforms_threads"]
    log:
        "logs/quantify/run_oarfish/{data_type}/{sample}.out",
    benchmark:
        repeat(
            "benchmarks/quantify/run_oarfish/{data_type}/{sample}.txt",
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


rule quantify_run_oarfish_novel:
    input:
        "results/align/run_minimap2_transcriptome_gencode_novel/{sample}/{sample}.aligned.bam",
    output:
        "results/quantify/run_oarfish_novel/{data_type}/{sample}/{sample}.quant",
    params:
        output_path="results/quantify/run_oarfish_novel/{data_type}/{sample}/{sample}",
        filter_group=config["oarfish_filter_group"],
        num_bootstraps=config["num_bootstraps"],
    threads: config["quantify_isoforms_threads"]
    log:
        "logs/quantify/run_oarfish_novel/{data_type}/{sample}.out",
    conda:
        "../envs/standalone/oarfish.yaml"
    shell:
        """
        conda list > {log};
        oarfish --threads {threads} \
                --filter-group {params.filter_group} \
                --model-coverage \
                --alignments {input} \
                --num-bootstraps {params.num_bootstraps} \
                --output {params.output_path} &>> {log}
        """


rule quantify_run_oarfish_bootstrapped:
    input:
        "results/align/run_minimap2_transcriptome_{data_type}/{sample}/{sample}.aligned.bam",
    output:
        "results/quantify/run_oarfish_bootstrapped/{data_type}/{sample}/{sample}.quant",
    params:
        output_path="results/quantify/run_oarfish_bootstrapped/{data_type}/{sample}/{sample}",
        filter_group=config["oarfish_filter_group"],
    threads: config["quantify_isoforms_threads"]
    log:
        "logs/quantify/run_oarfish_bootstrapped/{data_type}/{sample}.out",
    conda:
        "../envs/standalone/oarfish.yaml"
    shell:
        """
        conda list > {log};
        oarfish --threads {threads} \
                --filter-group {params.filter_group} \
                --model-coverage \
                --num-bootstraps 20 \
                --alignments {input} \
                --output {params.output_path} &>> {log}
        """


rule quantify_run_salmon:
    input:
        reads="results/align/run_minimap2_transcriptome_{data_type}/{sample}/{sample}.aligned.bam",
        transcriptome="results/prepare/extract_transcriptomes/{data_type}_transcriptome.fa",
    output:
        "results/quantify/run_salmon/{data_type}/{sample}/quant.sf",
    params:
        output_path="results/quantify/run_salmon/{data_type}/{sample}",
    threads: config["quantify_isoforms_threads"]
    log:
        "logs/quantify/run_salmon/{data_type}/{sample}.out",
    benchmark:
        repeat(
            "benchmarks/quantify/run_salmon/{data_type}/{sample}.txt",
            config["timing_repetitions"],
        )
    conda:
        "../envs/standalone/salmon.yaml"
    shell:
        """
        salmon --version > {log};
        salmon quant --ont -t {input.transcriptome} -a {input.reads} \
                -o {params.output_path} -p {threads} -l A \
                &>> {log}
        """


rule quantify_run_salmon_novel:
    input:
        reads="results/align/run_minimap2_transcriptome_gencode/{sample}/{sample}.aligned.bam",
        transcriptome="results/prepare/extract_novel_transcriptome/gencode_transcriptome.fa",
    output:
        "results/quantify/run_salmon_novel/gencode/{sample}/quant.sf",
    params:
        output_path="results/quantify/run_salmon_novel/gencode/{sample}",
    threads: config["quantify_isoforms_threads"]
    log:
        "logs/quantify/run_salmon_novel/gencode/{sample}.out",
    benchmark:
        repeat(
            "benchmarks/quantify/run_salmon_novel/gencode/{sample}.txt",
            config["timing_repetitions"],
        )
    conda:
        "../envs/standalone/salmon.yaml"
    shell:
        """
        salmon --version > {log};
        salmon quant --ont -t {input.transcriptome} -a {input.reads} \
                -o {params.output_path} -p {threads} -l A \
                &>> {log}
        """


rule quantify_run_bambu:
    input:
        reads="results/align/run_minimap2_{data_type}/{sample}/{sample}.aligned.sorted.bam",
        gencode_transcriptome="results/prepare/adjust_transcriptome_assembly_names/gencode.v45.primary_assembly.annotation.named.gtf",
        gencode_genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        sirv_transcriptome="results/prepare/adjust_sirv_names/sirv_set_four.gtf",
        sirv_genome="results/prepare/adjust_sirv_names/sirv_set_four.fa",
    output:
        "results/quantify/run_bambu/{data_type}/{sample}/counts_transcript.txt",
    threads: config["quantify_isoforms_threads"]
    params:
        outdir="results/quantify/run_bambu/{data_type}/{sample}",
        sirv="{data_type}",
        novel=False,
    log:
        "logs/quantify/run_bambu/{data_type}/{sample}.out",
    benchmark:
        repeat(
            "benchmarks/quantify/run_bambu/{data_type}/{sample}.txt",
            config["timing_repetitions"],
        )
    conda:
        "../envs/r/bambu.yaml"
    script:
        "../scripts/r/bambu_run_quantification.R"


rule quantify_run_bambu_novel:
    input:
        reads="results/align/run_minimap2_gencode/{sample}/{sample}.aligned.sorted.bam",
        gencode_genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        gencode_novel_transcriptome="results/discover/fix_bambu_gene_ids/transcriptome.gtf",
    output:
        "results/quantify/run_bambu_novel/gencode/{sample}/counts_transcript.txt",
    threads: config["quantify_isoforms_threads"]
    params:
        outdir="results/quantify/run_bambu_novel/gencode/{sample}",
        sirv="gencode",
        novel=True,
    log:
        "logs/quantify/run_bambu_novel/gencode/{sample}.out",
    benchmark:
        repeat(
            "benchmarks/quantify/run_bambu_novel/gencode/{sample}.txt",
            config["timing_repetitions"],
        )
    conda:
        "../envs/r/bambu.yaml"
    script:
        "../scripts/r/bambu_run_quantification.R"


rule quantify_run_isoquant:
    input:
        reads="results/align/run_minimap2_{data_type}/{sample}/{sample}.aligned.sorted.bam",
        gencode_genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        sirv_transcriptome="results/prepare/make_db_files/sirv_set_four.db",
        gencode_transcriptome="results/prepare/make_db_files/gencode.v45.primary_assembly.annotation.named.db",
        sirv_genome="results/prepare/adjust_sirv_names/sirv_set_four.fa",
    output:
        "results/quantify/run_isoquant/{data_type}/{sample}/OUT/OUT.transcript_counts.tsv",
    threads: config["quantify_isoforms_threads"]
    params:
        outdir="results/quantify/run_isoquant/{data_type}/{sample}",
    log:
        "logs/quantify/run_isoquant/{data_type}/{sample}.out",
    benchmark:
        repeat(
            "benchmarks/quantify/run_isoquant/{data_type}/{sample}.txt",
            config["timing_repetitions"],
        )
    conda:
        "../envs/py/isoquant.yaml"
    shell:
        """
        isoquant.py --version > {log};
        if [ "{wildcards.data_type}" = "sirv" ]; then
            isoquant.py --reference {input.sirv_genome} \
                --genedb {input.sirv_transcriptome} \
                --bam {input.reads} \
                --polya_requirement never \
                --no_model_construction \
                --complete_genedb \
                --threads {threads} \
                --data_type pacbio_ccs -o {params.outdir} &>> {log}
        else
            isoquant.py --reference {input.gencode_genome} \
                --genedb {input.gencode_transcriptome} \
                --bam {input.reads} \
                --polya_requirement never \
                --no_model_construction \
                --complete_genedb \
                --threads {threads} \
                --data_type pacbio_ccs -o {params.outdir} &>> {log}
        fi
        """


rule quantify_run_isoquant_novel:
    input:
        reads="results/align/run_minimap2_gencode/{sample}/{sample}.aligned.sorted.bam",
        gencode_genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        sirv_transcriptome="results/prepare/make_db_files/sirv_set_four.db",
        gencode_transcriptome="results/prepare/make_db_files/gencode.v45.primary_assembly.annotation.named.db",
        gencode_novel_transcriptome="results/prepare/make_db_files/gencode.novel.db",
        sirv_genome="results/prepare/adjust_sirv_names/sirv_set_four.fa",
    output:
        #directory("results/quantify/run_isoquant_novel/gencode/{sample}"),
        "results/quantify/run_isoquant_novel/gencode/{sample}/OUT/OUT.transcript_counts.tsv",
    threads: config["quantify_isoforms_threads"]
    params:
        outdir="results/quantify/run_isoquant_novel/gencode/{sample}",
    log:
        "logs/quantify/run_isoquant_novel/gencode/{sample}.out",
    benchmark:
        repeat(
            "benchmarks/quantify/run_isoquant_novel/gencode/{sample}.txt",
            config["timing_repetitions"],
        )
    conda:
        "../envs/py/isoquant.yaml"
    shell:
        """
        isoquant.py --reference {input.gencode_genome} \
            --genedb {input.gencode_novel_transcriptome} \
            --bam {input.reads} \
            --polya_requirement never \
            --no_model_construction \
            --complete_genedb \
            --threads {threads} \
            --data_type pacbio_ccs -o {params.outdir} &>> {log}
        """


rule quantify_run_kallisto_long:
    input:
        reads="results/prepare/convert_mapped_bams_to_fastq/{sample}/{data_type}.reads.fastq.gz",
        index=f"results/prepare/create_lr_kallisto_index/{{data_type}}_k-{config['lr_kallisto_index_k']}.idx",
        kallisto_binary="results/prepare/compile_lr_kallisto",
        gencode_transcriptome_gmap_headered="results/prepare/standardize_gtf_files/gencode_map.txt",
        sirv_four_transcriptome_gmap_headered="results/prepare/standardize_gtf_files/sirv_set_four_map.txt",
    output:
        "results/quantify/run_kallisto_long/{data_type}/{sample}/abundance_1.tsv",
    threads: config["quantify_isoforms_threads"]
    params:
        outdir="results/quantify/run_kallisto_long/{data_type}/{sample}",
        bus_threshold=config["kallisto_bus_threshold"],
        num_bootstraps=config["num_bootstraps"],
        seed=config["seed"],
    log:
        "logs/quantify/run_kallisto_long/{data_type}/{sample}.out",
    benchmark:
        repeat(
            "benchmarks/quantify/run_kallisto_long/{data_type}/{sample}.txt",
            config["timing_repetitions"],
        )
    conda:
        "../envs/standalone/kallisto.yaml"
    shell:
        """
        touch {log};
        if [ "{wildcards.data_type}" = "sirv" ]; then
            {input.kallisto_binary}/kallisto/build/src/kallisto  bus -t {threads} --long --threshold \
                {params.bus_threshold} -x bulk -i {input.index} \
                -o {params.outdir} {input.reads} &>> {log};
            bustools sort -t {threads} {params.outdir}/output.bus \
                -o {params.outdir}/sorted.bus &>> {log}; 
            bustools count {params.outdir}/sorted.bus \
            -t {params.outdir}/transcripts.txt \
            -e {params.outdir}/matrix.ec \
            -g {input.sirv_four_transcriptome_gmap_headered} \
            -o {params.outdir}/count --cm -m \
            &>> {log};
            {input.kallisto_binary}/kallisto/build/src/kallisto  quant-tcc -t {threads} \
                --long -P PacBio -f {params.outdir}/flens.txt \
                {params.outdir}/count.mtx -i {input.index} \
                -e {params.outdir}/count.ec.txt \
                -o {params.outdir} \
                -b {params.num_bootstraps} \
                --seed {params.seed} \
                --matrix-to-files \
                &>> {log}
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
                -o {params.outdir} \
                -b {params.num_bootstraps} \
                --seed {params.seed} \
                --matrix-to-files \
                &>> {log}
        fi
        """


rule quantify_run_kallisto_long_novel:
    input:
        reads="results/prepare/convert_mapped_bams_to_fastq/{sample}/gencode.reads.fastq.gz",
        index=f"results/prepare/create_lr_kallisto_index_novel/gencode_novel_k-{config['lr_kallisto_index_k']}.idx",
        kallisto_binary="results/prepare/compile_lr_kallisto",
        gencode_novel_transcriptome_gmap_headered="results/prepare/standardize_gtf_files_novel/gencode_novel_map.txt",
    output:
        #directory("results/quantify/run_kallisto_long_novel/gencode/{sample}"),
        "results/quantify/run_kallisto_long_novel/gencode/{sample}/abundance_1.tsv",
    threads: config["quantify_isoforms_threads"]
    params:
        outdir="results/quantify/run_kallisto_long_novel/gencode/{sample}",
        bus_threshold=config["kallisto_bus_threshold"],
        num_bootstraps=config["num_bootstraps"],
        seed=config["seed"],
    log:
        "logs/quantify/run_kallisto_long_novel/gencode/{sample}.out",
    conda:
        "../envs/standalone/kallisto.yaml"
    shell:
        """
        touch {log};
        {input.kallisto_binary}/kallisto/build/src/kallisto  bus -t {threads} --long --threshold \
            {params.bus_threshold} -x bulk -i {input.index} \
            -o {params.outdir} {input.reads} &>> {log};
        bustools sort -t {threads} {params.outdir}/output.bus \
            -o {params.outdir}/sorted.bus &>> {log}; 
        bustools count {params.outdir}/sorted.bus \
        -t {params.outdir}/transcripts.txt \
        -e {params.outdir}/matrix.ec \
        -g {input.gencode_novel_transcriptome_gmap_headered} \
        -o {params.outdir}/count --cm -m \
        &>> {log};
        {input.kallisto_binary}/kallisto/build/src/kallisto  quant-tcc -t {threads} \
            --long -P PacBio -f {params.outdir}/flens.txt \
            {params.outdir}/count.mtx -i {input.index} \
            -e {params.outdir}/count.ec.txt \
            -o {params.outdir} \
            -b {params.num_bootstraps} \
            --seed {params.seed} \
            --matrix-to-files \
            &>> {log}

        """


rule quantify_run_salmon_illumina:
    input:
        index="results/prepare/create_salmon_index/{data_type}",
        first_reads="results/align/convert_mapped_bams_to_fastq/{sample}/{data_type}/{sample}.reads_1.fastq.gz",
        second_reads="results/align/convert_mapped_bams_to_fastq/{sample}/{data_type}/{sample}.reads_2.fastq.gz",
    output:
        "results/quantify/run_salmon_illumina/{data_type}/{sample}/quant.sf",
    params:
        output_name="results/quantify/run_salmon_illumina/{data_type}/{sample}",
        fragment_length_mean=config["salmon_fragment_length_mean"],
        fragment_length_sd=config["salmon_fragment_length_sd"],
        num_bootstraps=config["num_bootstraps"],
    log:
        "logs/quantify/run_salmon_illumina/{data_type}/{sample}.out",
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
            --fldSD {params.fragment_length_sd} &>> {log} \
            --numBootstraps {params.num_bootstraps} &>> {log}
        """


rule quantify_run_salmon_illumina_novel:
    input:
        index="results/prepare/create_salmon_index_novel/gencode",
        first_reads="results/align/convert_mapped_bams_to_fastq/{sample}/gencode/{sample}.reads_1.fastq.gz",
        second_reads="results/align/convert_mapped_bams_to_fastq/{sample}/gencode/{sample}.reads_2.fastq.gz",
    output:
        "results/quantify/run_salmon_illumina_novel/gencode/{sample}/quant.sf",
    params:
        output_name="results/quantify/run_salmon_illumina_novel/gencode/{sample}",
        fragment_length_mean=config["salmon_fragment_length_mean"],
        fragment_length_sd=config["salmon_fragment_length_sd"],
        num_bootstraps=config["num_bootstraps"],
    log:
        "logs/quantify/run_salmon_illumina_novel/gencode/{sample}.out",
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
            --fldSD {params.fragment_length_sd} &>> {log} \
            --numBootstraps {params.num_bootstraps} &>> {log}
        """
