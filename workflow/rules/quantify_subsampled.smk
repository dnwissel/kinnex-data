configfile: "config/config.yaml"


singularity: f"docker://condaforge/mambaforge:{config['mambaforge_version']}"


rule quantify_subsampled_run_oarfish:
    input:
        "results/downsample/run_minimap2_transcriptome_{type}/{subsample_number}/{read_number}/{sample}/{sample}.aligned.bam",
    output:
        "results/quantify_subsampled/run_oarfish/{type}/{subsample_number}/{read_number}/{sample}/{sample}.quant",
    params:
        output_path="results/quantify_subsampled/run_oarfish/{type}/{subsample_number}/{read_number}/{sample}/{sample}",
        filter_group=config["oarfish_filter_group"],
    benchmark:
        "benchmarks/quantify_subsampled/run_oarfish/{type}/{subsample_number}/{read_number}/{sample}.txt"
    log:
        "logs/quantify_subsampled/run_oarfish/{type}/{subsample_number}/{read_number}/{sample}.out",
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


rule quantify_subsampled_run_oarfish_bootstrap:
    input:
        "results/downsample/run_minimap2_transcriptome_{type}/{subsample_number}/{read_number}/{sample}/{sample}.aligned.bam",
    output:
        "results/quantify_subsampled/run_oarfish_bootstrap/{type}/{subsample_number}/{read_number}/{sample}/{sample}.quant",
    params:
        output_path="results/quantify_subsampled/run_oarfish_bootstrap/{type}/{subsample_number}/{read_number}/{sample}/{sample}",
        filter_group=config["oarfish_filter_group"],
        num_bootstraps=config["num_bootstraps"],
    benchmark:
        "benchmarks/quantify_subsampled/run_oarfish_bootstrap/{type}/{subsample_number}/{read_number}/{sample}.txt"
    log:
        "logs/quantify_subsampled/run_oarfish_bootstrap/{type}/{subsample_number}/{read_number}/{sample}.out",
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
                --num-bootstraps {params.num_bootstraps} \
                --output {params.output_path} &>> {log}
        """


rule quantify_subsampled_run_salmon_illumina:
    input:
        index="results/prepare/create_salmon_index/{type}",
        first_reads="results/downsample/convert_mapped_bams_to_fastq/{subsample_number}_{read_number}_{type}/{sample}-r1.fastq.gz",
        second_reads="results/downsample/convert_mapped_bams_to_fastq/{subsample_number}_{read_number}_{type}/{sample}-r2.fastq.gz",
    output:
        "results/quantify_subsampled/run_salmon_illumina/{type}/{subsample_number}/{read_number}/{sample}/quant.sf",
    params:
        output_name="results/quantify_subsampled/run_salmon_illumina/{type}/{subsample_number}/{read_number}/{sample}",
        fragment_length_mean=config["salmon_fragment_length_mean"],
        fragment_length_sd=config["salmon_fragment_length_sd"],
    benchmark:
        "benchmarks/quantify_subsampled/run_salmon_illumina/{type}/{subsample_number}/{read_number}/{sample}.txt"
    log:
        "logs/quantify_subsampled/run_salmon_illumina/{type}/{subsample_number}/{read_number}/{sample}.out",
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


rule quantify_subsampled_run_salmon_illumina_bootstrap:
    input:
        index="results/prepare/create_salmon_index/{type}",
        first_reads="results/downsample/convert_mapped_bams_to_fastq/{subsample_number}_{read_number}_{type}/{sample}-r1.fastq.gz",
        second_reads="results/downsample/convert_mapped_bams_to_fastq/{subsample_number}_{read_number}_{type}/{sample}-r2.fastq.gz",
    output:
        "results/quantify_subsampled/run_salmon_illumina_bootstrap/{type}/{subsample_number}/{read_number}/{sample}/quant.sf",
    params:
        output_name="results/quantify_subsampled/run_salmon_illumina_bootstrap/{type}/{subsample_number}/{read_number}/{sample}",
        fragment_length_mean=config["salmon_fragment_length_mean"],
        fragment_length_sd=config["salmon_fragment_length_sd"],
        num_bootstraps=config["num_bootstraps"],
    log:
        "logs/quantify_subsampled/run_salmon_illumina_bootstrap/{type}/{subsample_number}/{read_number}/{sample}.out",
    threads: config["quantify_isoforms_threads"]
    conda:
        "../envs/standalone/salmon.yaml"
    shell:
        """
        salmon --version > {log};
        salmon quant -i {input.index} -l A \
            -1 {input.first_reads} -2 {input.second_reads} \
            --validateMappings -o {params.output_name} -p {threads} \
            --numBootstraps {params.num_bootstraps} \
            --seqBias --gcBias --fldMean {params.fragment_length_mean} \
            --fldSD {params.fragment_length_sd} &>> {log}
        """


rule quantify_subsampled_run_kallisto_long_bootstrap:
    input:
        reads="results/downsample/convert_mapped_bams_to_fastq_long_reads/{subsample_number}_{read_number}_{type}/{sample}.fastq.gz",
        index=f"results/prepare/create_lr_kallisto_index/{{type}}_k-{config['lr_kallisto_index_k']}.idx",
        kallisto_binary="results/prepare/compile_lr_kallisto",
        gencode_transcriptome_gmap_headered="results/prepare/standardize_gtf_files/gencode_map.txt",
        sirv_four_transcriptome_gmap_headered="results/prepare/standardize_gtf_files/sirv_set_four_map.txt",
    output:
        #directory(
        #    "results/quantify_subsampled/run_kallisto_long/{type}/{subsample_number}/{read_number}/{sample}"
        #),
        "results/quantify_subsampled/run_kallisto_long_bootstrap/{type}/{subsample_number}/{read_number}/{sample}/abundance_1.tsv",
        #"results/quantify/run_kallisto_long_novel/gencode/{sample}/abundance.tsv",
    threads: config["quantify_isoforms_threads"]
    params:
        outdir="results/quantify_subsampled/run_kallisto_long_bootstrap/{type}/{subsample_number}/{read_number}/{sample}",
        bus_threshold=config["kallisto_bus_threshold"],
        num_bootstraps=config["num_bootstraps"],
        seed=config["seed"],
    log:
        "logs/quantify_subsampled/run_kallisto_long_bootstrap/{type}/{subsample_number}/{read_number}/{sample}.out",
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


rule quantify_subsampled_run_salmon:
    input:
        reads="results/downsample/run_minimap2_transcriptome_{type}/{subsample_number}/{read_number}/{sample}/{sample}.aligned.bam",
        transcriptome="results/prepare/extract_transcriptomes/{type}_transcriptome.fa",
    output:
        "results/quantify_subsampled/run_salmon/{type}/{subsample_number}/{read_number}/{sample}/quant.sf",
    params:
        output_path="results/quantify_subsampled/run_salmon/{type}/{subsample_number}/{read_number}/{sample}",
    benchmark:
        "benchmarks/quantify_subsampled/run_salmon/{type}/{subsample_number}/{read_number}/{sample}.txt"
    log:
        "logs/quantify_subsampled/run_salmon/{type}/{subsample_number}/{read_number}/{sample}.out",
    threads: config["quantify_isoforms_threads"]
    conda:
        "../envs/standalone/salmon.yaml"
    shell:
        """
        salmon --version > {log};
        salmon quant --ont -t {input.transcriptome} -a {input.reads} \
                -o {params.output_path} -p {threads} -l A &>> {log}
        """


rule quantify_subsampled_run_bambu:
    input:
        reads="results/downsample/run_minimap2_{type}/{subsample_number}/{read_number}/{sample}/{sample}.aligned.sorted.bam",
        gencode_transcriptome="results/prepare/adjust_transcriptome_assembly_names/gencode.v45.primary_assembly.annotation.named.gtf",
        gencode_genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        sirv_transcriptome="results/prepare/adjust_sirv_names/sirv_set_four.gtf",
        sirv_genome="results/prepare/adjust_sirv_names/sirv_set_four.fa",
    output:
        "results/quantify_subsampled/run_bambu/{type}/{subsample_number}/{read_number}/{sample}/counts_transcript.txt",
    threads: config["quantify_isoforms_threads"]
    params:
        outdir="results/quantify_subsampled/run_bambu/{type}/{subsample_number}/{read_number}/{sample}",
        sirv="{type}",
        novel=False,
    benchmark:
        "benchmarks/quantify_subsampled/run_bambu/{type}/{subsample_number}/{read_number}/{sample}.txt"
    log:
        "logs/quantify_subsampled/run_bambu/{type}/{subsample_number}/{read_number}/{sample}.out",
    conda:
        "../envs/r/bambu.yaml"
    script:
        "../scripts/r/bambu_run_quantification.R"


rule quantify_subsampled_run_isoquant:
    input:
        reads="results/downsample/run_minimap2_{type}/{subsample_number}/{read_number}/{sample}/{sample}.aligned.sorted.bam",
        gencode_genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        sirv_transcriptome="results/prepare/make_db_files/sirv_set_four.db",
        gencode_transcriptome="results/prepare/make_db_files/gencode.v45.primary_assembly.annotation.named.db",
        sirv_genome="results/prepare/adjust_sirv_names/sirv_set_four.fa",
    output:
        #directory(
        #    "results/quantify_subsampled/run_isoquant/{type}/{subsample_number}/{read_number}/{sample}"
        #),
        "results/quantify_subsampled/run_isoquant/{type}/{subsample_number}/{read_number}/{sample}/OUT/OUT.transcript_counts.tsv",
    threads: config["quantify_isoforms_threads"]
    params:
        outdir="results/quantify_subsampled/run_isoquant/{type}/{subsample_number}/{read_number}/{sample}",
    benchmark:
        "benchmarks/quantify_subsampled/run_isoquant/{type}/{subsample_number}/{read_number}/{sample}.txt"
    log:
        "logs/quantify_subsampled/run_isoquant/{type}/{subsample_number}/{read_number}/{sample}.out",
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


rule quantify_subsampled_run_kallisto_long:
    input:
        reads="results/downsample/convert_mapped_bams_to_fastq_long_reads/{subsample_number}_{read_number}_{type}/{sample}.fastq.gz",
        index=f"results/prepare/create_lr_kallisto_index/{{type}}_k-{config['lr_kallisto_index_k']}.idx",
        kallisto_binary="results/prepare/compile_lr_kallisto",
        gencode_transcriptome_gmap_headered="results/prepare/standardize_gtf_files/gencode_map.txt",
        sirv_four_transcriptome_gmap_headered="results/prepare/standardize_gtf_files/sirv_set_four_map.txt",
    output:
        "results/quantify_subsampled/run_kallisto_long/{type}/{subsample_number}/{read_number}/{sample}/abundance_1.tsv",
    threads: config["quantify_isoforms_threads"]
    params:
        outdir="results/quantify_subsampled/run_kallisto_long/{type}/{subsample_number}/{read_number}/{sample}",
        bus_threshold=config["kallisto_bus_threshold"],
    benchmark:
        "benchmarks/quantify_subsampled/run_kallisto_long/{type}/{subsample_number}/{read_number}/{sample}.txt"
    log:
        "logs/quantify_subsampled/run_kallisto_long/{type}/{subsample_number}/{read_number}/{sample}.out",
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
            -g {input.sirv_four_transcriptome_gmap_headered} \
            -o {params.outdir}/count --cm -m \
            &>> {log};
            {input.kallisto_binary}/kallisto/build/src/kallisto  quant-tcc -t {threads} \
                --long -P PacBio -f {params.outdir}/flens.txt \
                {params.outdir}/count.mtx -i {input.index} \
                -e {params.outdir}/count.ec.txt \
                -o {params.outdir} \
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
                --matrix-to-files \
                &>> {log}
        fi
        """
