configfile: "config/config.yaml"


rule prepare_download_genome:
    output:
        "results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    params:
        url=config["genome_url"],
        download_path="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
    log:
        "logs/prepare/download_genome/out.log",
    conda:
        "../envs/standalone/curl.yaml"
    shell:
        """
        wget -O {params.download_path} {params.url} &> {log};
        gunzip {params.download_path}
        """


rule prepare_download_transcriptome:
    output:
        "results/prepare/download_transcriptome/gencode.v45.primary_assembly.annotation.gtf.gz",
    params:
        url=config["transcriptome_url"],
    log:
        "logs/prepare/download_transcriptome/out.log",
    shell:
        """
        wget -O {output} {params.url} &> {log}
        """


rule prepare_download_sqanti_polya_motif_list:
    output:
        "results/prepare/download_sqanti_polya_motif_list/motifs.txt",
    params:
        url=config["polya_motif_url"],
    log:
        "logs/prepare/download_sqanti_polya_motif_list/out.log",
    shell:
        """
        wget -O {output} {params.url} &> {log}
        """


rule prepare_download_sirvome:
    output:
        sirv_set_four=directory(
            "results/prepare/download_sirvome/SIRV_Set4_Norm_Sequences_20210507"
        ),
    params:
        set_four_url=config["sirv_set_four_url"],
        set_four_download_object="results/prepare/download_sirvome/set_four.zip",
        download_path="results/prepare/download_sirvome/",
    log:
        "logs/prepare/download_sirvome/out.log",
    conda:
        "../envs/standalone/curl.yaml"
    shell:
        """
        curl -L {params.set_four_url} \
            --output {params.set_four_download_object} &> {log};
        unzip -o {params.set_four_download_object} \
            -d {params.download_path} >> {log}
        """


rule prepare_adjust_transcriptome_assembly_names:
    input:
        genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        transcriptome="results/prepare/download_transcriptome/gencode.v45.primary_assembly.annotation.gtf.gz",
    output:
        "results/prepare/adjust_transcriptome_assembly_names/gencode.v45.primary_assembly.annotation.named.gtf",
    log:
        "logs/prepare/adjust_transcriptome_assembly_names/out.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/prepare/adjust_transcriptome_assembly_names.R"


rule prepare_adjust_sirv_names:
    input:
        "results/prepare/download_sirvome/SIRV_Set4_Norm_Sequences_20210507",
    output:
        output_path_transcriptome="results/prepare/adjust_sirv_names/sirv_set_four.gtf",
        output_path_genome="results/prepare/adjust_sirv_names/sirv_set_four.fa",
    params:
        transcriptome_name=config["sirv_transcriptome_name"],
        genome_name=config["sirv_genome_name"],
    log:
        "logs/prepare/adjust_sirv_names/out.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/prepare/adjust_sirv_names.R"


rule prepare_shuffle_fastq:
    input:
        "results/prepare/convert_mapped_bams_to_fastq/{sample}/{type}.reads.fastq.gz",
    output:
        "results/prepare/shuffle_fastq/{sample}/{type}.reads.fastq.gz",
    params:
        ziplevel=9,
        mem_gb=12,
    log:
        "logs/prepare/shuffle_fastq/{type}/{sample}/out.log",
    threads: 12
    conda:
        "../envs/standalone/bbmap.yaml"
    shell:
        """
        shuffle.sh in={input} out={output} ziplevel={params.ziplevel} -Xmx{params.mem_gb}g
        """


rule prepare_standardize_gtf_files:
    input:
        sirv_four_transcriptome="results/prepare/adjust_sirv_names/sirv_set_four.gtf",
        gencode_transcriptome="results/prepare/adjust_transcriptome_assembly_names/gencode.v45.primary_assembly.annotation.named.gtf",
    output:
        sirv_four_transcriptome="results/prepare/standardize_gtf_files/sirv_set_four.gtf",
        sirv_gff="results/prepare/standardize_gtf_files/sirv_set_four.gff",
        gencode_transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        gencode_gff="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gff",
        gencode_transcriptome_gmap="results/prepare/standardize_gtf_files/gencode_map.txt",
        gencode_transcriptome_gmap_headered="results/prepare/standardize_gtf_files/gencode_map_headered.txt",
        sirv_four_transcriptome_gmap="results/prepare/standardize_gtf_files/sirv_set_four_map.txt",
        sirv_four_transcriptome_gmap_headered="results/prepare/standardize_gtf_files/sirv_set_four_map_headered.txt",
    log:
        "logs/prepare/standardize_gtf_files/out.log",
    conda:
        "../envs/standalone/gffread.yaml"
    shell:
        """
        gffread --version > {log};
        gffread -E {input.sirv_four_transcriptome} \
            -T -o {output.sirv_four_transcriptome} &>> {log};
        gffread -E {input.gencode_transcriptome} \
            -T -o {output.gencode_transcriptome} &>> {log};
        gffread {output.sirv_four_transcriptome} \
            --table transcript_id,gene_id \
            > {output.sirv_four_transcriptome_gmap} 2>> {log};
        gffread {output.gencode_transcriptome} \
            --table transcript_id,gene_id \
            > {output.gencode_transcriptome_gmap} 2>> {log};
        echo -e "transcript\tgene" | \
            cat - {output.gencode_transcriptome_gmap} \
                > {output.gencode_transcriptome_gmap_headered} 2>> {log};
        echo -e "transcript\tgene" | \
            cat - {output.sirv_four_transcriptome_gmap} \
                > {output.sirv_four_transcriptome_gmap_headered} 2>> {log};
        gffread -o {output.gencode_gff} {output.gencode_transcriptome} &>> {log};
        gffread -o {output.sirv_gff} {output.sirv_four_transcriptome} &>> {log}
        """


rule prepare_make_db_files_sirv:
    input:
        input_gtf="results/prepare/standardize_gtf_files/sirv_set_four.gtf",
    output:
        output_db="results/prepare/make_db_files/sirv_set_four.db",
    params:
        checklines=config["gffutils_checklines"],
    log:
        "logs/prepare/make_db_files/sirv.log",
    conda:
        "../envs/standalone/gffutils.yaml"
    script:
        "../scripts/py/create_db_files.py"


rule prepare_make_db_files_gencode:
    input:
        input_gtf="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
    output:
        output_db="results/prepare/make_db_files/gencode.v45.primary_assembly.annotation.named.db",
    params:
        checklines=config["gffutils_checklines"],
    log:
        "logs/prepare/make_db_files/gencode.log",
    conda:
        "../envs/standalone/gffutils.yaml"
    script:
        "../scripts/py/create_db_files.py"


rule prepare_make_db_files_gencode_novel:
    input:
        input_gtf="results/prepare/standardize_gtf_files_novel/gencode_novel.gtf",
    output:
        output_db="results/prepare/make_db_files/gencode.novel.db",
    params:
        checklines=config["gffutils_checklines"],
    log:
        "logs/prepare/make_db_files/novel.log",
    conda:
        "../envs/standalone/gffutils.yaml"
    script:
        "../scripts/py/create_db_files.py"


rule prepare_concatenate_genomes:
    input:
        sirv_four_genome="results/prepare/adjust_sirv_names/sirv_set_four.fa",
        gencode_genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    output:
        overall_genome="results/prepare/concatenate_genomes/genome.fa",
    conda:
        "../envs/standalone/samtools.yaml"
    log:
        "logs/prepare/concatenate_genomes/out.log",
    shell:
        """
        conda list &> {log};
        cat {input.gencode_genome} {input.sirv_four_genome} > \
            {output.overall_genome} 2> {log};
        samtools faidx {output.overall_genome} &>> {log}
        """


rule prepare_concatenate_transcriptomes:
    input:
        sirv_four_transcriptome="results/prepare/standardize_gtf_files/sirv_set_four.gtf",
        gencode_transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
    output:
        overall_transcriptome="results/prepare/concatenate_transcriptomes/transcriptome.gtf",
    log:
        "logs/prepare/concatenate_transcriptomes/out.log",
    shell:
        """
        cat {input.gencode_transcriptome} {input.sirv_four_transcriptome} \
            > {output.overall_transcriptome} 2> {log}
        """


rule prepare_convert_ubam_to_fastqz:
    input:
        "data/seq/kinnex/raw/{sample}/flnc.bam",
    output:
        "results/prepare/convert_ubam_to_fastqz/{sample}/read.fastq.gz",
    params:
        compression_level=config["compression_level"],
    threads: config["stall_io_threads"]
    log:
        "logs/prepare/convert_ubam_to_fastqz/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    shell:
        """
        samtools --version > {log};
        samtools fastq -@ {threads} -0 {output} \
            -c {params.compression_level} {input} &>> {log};
        rm {input} &>> {log}
        """


rule prepare_run_minimap2:
    input:
        reads="results/prepare/convert_ubam_to_fastqz/{sample}/read.fastq.gz",
        genome="results/prepare/concatenate_genomes/genome.fa",
    output:
        "results/prepare/run_minimap2/{sample}/aligned.sorted.bam",
    params:
        align_sort_bam_memory_gb=config["align_sort_bam_memory_gb"],
        align_sort_bam_threads=config["align_sort_bam_threads"],
    threads: config["align_map_bam_threads"]
    log:
        "logs/prepare/run_minimap2/{sample}.log",
    conda:
        "../envs/standalone/minimap2.yaml"
    shell:
        """
        minimap2 --version > {log};
        minimap2 -ax splice:hq -t {threads} -uf {input.genome} \
            {input.reads} 2>> {log} | \
            samtools sort -@ {params.align_sort_bam_threads} \
            -m{params.align_sort_bam_memory_gb}g -o {output} - &>> {log};
        sleep 61s &>> {log};
        samtools index {output} &>> {log};
        rm {input.reads} &>> {log}
        """


rule prepare_filter_genome_mappings:
    input:
        "results/prepare/run_minimap2/{sample}/aligned.sorted.bam",
    output:
        gencode_reads="results/prepare/filter_genome_mappings/{sample}/{sample}.aligned.gencode.sorted.bam",
        gencode_ix="results/prepare/filter_genome_mappings/{sample}/{sample}.aligned.gencode.sorted.bam.bai",
        sirv_reads="results/prepare/filter_genome_mappings/{sample}/{sample}.aligned.sirv.sorted.bam",
        sirv_ix="results/prepare/filter_genome_mappings/{sample}/{sample}.aligned.sirv.sorted.bam.bai",
    params:
        lines_to_cut=1,
        gencode_chromosomes="chr",
    threads: config["stall_io_threads"]
    log:
        "logs/prepare/filter_genome_mappings/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    shell:
        """
        conda list &> {log};
        samtools idxstats {input} 2>> {log} | \
            cut -f {params.lines_to_cut} 2>> {log} | \
            grep {params.gencode_chromosomes}  2>> {log} | \
            xargs samtools view -b {input} > {output.gencode_reads} 2>> {log};
        samtools idxstats {input} 2>> {log} | \
            cut -f {params.lines_to_cut} 2>> {log} | \
            grep -v {params.gencode_chromosomes} 2>> {log} | \
            head -n -{params.lines_to_cut} 2>> {log} | \
            xargs samtools view -b {input} > {output.sirv_reads} 2>> {log};
        sleep 61s &>> {log};
        samtools index {output.gencode_reads} &>> {log};
        samtools index {output.sirv_reads} &>> {log};
        rm {input} &>> {log}
        """


rule prepare_convert_mapped_bams_to_fastq:
    input:
        "results/prepare/filter_genome_mappings/{sample}/{sample}.aligned.{type}.sorted.bam",
    output:
        "results/prepare/convert_mapped_bams_to_fastq/{sample}/{type}.reads.fastq.gz",
    params:
        compression_level=config["compression_level"],
    threads: config["stall_io_threads"]
    log:
        "logs/prepare/filter_genome_mappings/{type}/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    shell:
        """
        samtools --version > {log};
        samtools fastq -@ {threads} -0 {output} \
            -c {params.compression_level} {input} &>> {log};
        rm {input} &>> {log}
        """


rule prepare_convert_mapped_bams_to_fastq_random:
    input:
        "results/prepare/filter_genome_mappings/{sample}/{sample}.aligned.{type}.sorted.bam",
    output:
        "results/prepare/convert_mapped_bams_to_fastq_random/{sample}/{type}.reads.fastq.gz",
    params:
        compression_level=config["compression_level"],
    threads: config["stall_io_threads"]
    log:
        "logs/prepare/convert_mapped_bams_to_fastq_random/{type}/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    shell:
        """
        samtools --version > {log};
        samtools sort -n -m{params.memory}g -@ {threads} {input} 2>> {log} | samtools fastq -@ {threads} -0 {output} \
            -c {params.compression_level} - &>> {log};
        rm {input} &>> {log}
        """


rule prepare_convert_gtfs_to_beds_gencode:
    input:
        gencode_gtf="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        sirv_gtf="results/prepare/standardize_gtf_files/sirv_set_four.gtf",
    output:
        gencode_bed="results/prepare/convert_gtfs_to_beds/gencode.bed",
        sirv_bed="results/prepare/convert_gtfs_to_beds/sirv.bed",
    log:
        "logs/prepare/convert_gtfs_to_beds/out.log",
    conda:
        "../envs/standalone/minimap2.yaml"
    shell:
        """
        paftools.js gff2bed {input.gencode_gtf} > {output.gencode_bed} 2>> {log};
        paftools.js gff2bed {input.sirv_gtf} > {output.sirv_bed} 2>> {log};
        """


rule prepare_extract_transcriptomes:
    input:
        gencode_genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        sirv_transcriptome="results/prepare/standardize_gtf_files/sirv_set_four.gtf",
        gencode_transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        joint_genome="results/prepare/concatenate_genomes/genome.fa",
        joint_transcriptome="results/prepare/concatenate_transcriptomes/transcriptome.gtf",
        sirv_genome="results/prepare/adjust_sirv_names/sirv_set_four.fa",
    output:
        gencode_transcriptome="results/prepare/extract_transcriptomes/gencode_transcriptome.fa",
        joint_transcriptome="results/prepare/extract_transcriptomes/joint_transcriptome.fa",
        sirv_transcriptome="results/prepare/extract_transcriptomes/sirv_transcriptome.fa",
    log:
        "logs/prepare/extract_transcriptomes/out.log",
    conda:
        "../envs/standalone/gffread.yaml"
    shell:
        """
        gffread -w {output.gencode_transcriptome} \
            -g {input.gencode_genome} \
            {input.gencode_transcriptome} &> {log};
        gffread -w {output.sirv_transcriptome} \
            -g {input.sirv_genome} \
            {input.sirv_transcriptome} &>> {log};
        gffread -w {output.joint_transcriptome} \
            -g {input.joint_genome} \
            {input.joint_transcriptome} &>> {log}
        """


rule prepare_standardize_gtf_files_novel:
    input:
        "results/discover/fix_bambu_gene_ids/transcriptome.gtf",
    output:
        transcriptome="results/prepare/standardize_gtf_files_novel/gencode_novel.gtf",
        gff="results/prepare/standardize_gtf_files_novel/gencode_novel.gff",
        transcriptome_gmap="results/prepare/standardize_gtf_files_novel/gencode_novel_map.txt",
        transcriptome_gmap_headered="results/prepare/standardize_gtf_files_novel/gencode_novel_map_headered.txt",
    log:
        "logs/prepare/standardize_gtf_files_novel/out.log",
    conda:
        "../envs/standalone/gffread.yaml"
    shell:
        """
        gffread --version > {log};
        gffread -E {input} \
            -T -o {output.transcriptome} &>> {log};
        gffread {output.transcriptome} \
            --table transcript_id,gene_id \
            > {output.transcriptome_gmap} 2>> {log};
        echo -e "transcript\tgene" | \
            cat - {output.transcriptome_gmap} \
                > {output.transcriptome_gmap_headered} 2>> {log};
        gffread -o {output.gff} {output.transcriptome} &>> {log}
        """


rule prepare_extract_novel_transcriptome:
    input:
        gencode_genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        gencode_novel_transcriptome="results/discover/fix_bambu_gene_ids/transcriptome.gtf",
    output:
        gencode_novel_transcriptome="results/prepare/extract_novel_transcriptome/gencode_transcriptome.fa",
    log:
        "logs/prepare/extract_novel_transcriptome/out.log",
    conda:
        "../envs/standalone/gffread.yaml"
    shell:
        """
        conda list &> {log};
        gffread -w {output.gencode_novel_transcriptome} \
            -g {input.gencode_genome} \
            {input.gencode_novel_transcriptome} &> {log}
        """


rule prepare_run_trim_galore:
    input:
        reads_first="data/seq/kinnex/raw/illumina/FASTQ/{sample}-r1.fastq.gz",
        reads_second="data/seq/kinnex/raw/illumina/FASTQ/{sample}-r2.fastq.gz",
    output:
        reads_first="results/prepare/run_trim_galore/{sample}-r1_val_1.fq.gz",
        reads_second="results/prepare/run_trim_galore/{sample}-r2_val_2.fq.gz",
    params:
        min_quality=config["illumina_trim_min_quality"],
        min_length=config["illumina_trim_min_length"],
        outdir="results/prepare/run_trim_galore",
    threads: config["stall_io_threads"]
    log:
        "logs/prepare/run_trim_galore/{sample}.log",
    conda:
        "../envs/standalone/trim_galore.yaml"
    shell:
        """
        trim_galore --version > {log};
        trim_galore -q {params.min_quality} --phred33 \
            --length {params.min_length} -o {params.outdir} \
            --paired {input.reads_first} {input.reads_second} &>> {log}
        """


rule prepare_index_star:
    input:
        genome="results/prepare/concatenate_genomes/genome.fa",
        transcriptome="results/prepare/concatenate_transcriptomes/transcriptome.gtf",
    output:
        directory("results/prepare/index_star"),
    params:
        read_length=config["illumina_read_length"],
    threads: config["align_map_bam_threads"]
    log:
        "logs/prepare/index_star/out.log",
    conda:
        "../envs/standalone/star.yaml"
    shell:
        """
        STAR --version > {log};
        STAR --runMode genomeGenerate --runThreadN {threads} \
            --genomeDir {output} --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.transcriptome} \
            --sjdbOverhang {params.read_length} &>> {log}
        """


rule prepare_create_decoys_salmon_index:
    input:
        sirv_genome="results/prepare/adjust_sirv_names/sirv_set_four.fa",
        gencode_genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        gencode_transcriptome="results/prepare/extract_transcriptomes/gencode_transcriptome.fa",
        sirv_transcriptome="results/prepare/extract_transcriptomes/sirv_transcriptome.fa",
    output:
        sirv_decoys="results/prepare/create_decoys_salmon_index/sirv_decoys.txt",
        gencode_decoys="results/prepare/create_decoys_salmon_index/gencode_decoys.txt",
        gencode_decoyed_transcriptome="results/prepare/create_decoys_salmon_index/gencode_decoyed_sequence.fa",
        sirv_decoyed_transcriptome="results/prepare/create_decoys_salmon_index/sirv_decoyed_sequence.fa",
    params:
        outdir="results/prepare/create_decoys_salmon_index",
    log:
        "logs/prepare/create_decoys_salmon_index/out.log",
    shell:
        """
        mkdir -p {params.outdir};
        grep "^>" {input.sirv_genome} | cut -d " " -f 1 \
            > {output.sirv_decoys} 2> {log}
        grep "^>" {input.gencode_genome} | cut -d " " -f 1 \
            > {output.gencode_decoys} 2>> {log};
        sed -i.bak -e 's/>//g' {output.sirv_decoys} 2>> {log};
        sed -i.bak -e 's/>//g' {output.gencode_decoys} 2>> {log};
        cat {input.gencode_transcriptome} {input.gencode_genome} \
            > {output.gencode_decoyed_transcriptome} 2>> {log};
        cat {input.sirv_transcriptome} {input.sirv_genome} \
            > {output.sirv_decoyed_transcriptome} 2>> {log}
        """


rule prepare_create_decoys_salmon_index_novel:
    input:
        gencode_genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        gencode_transcriptome="results/prepare/extract_novel_transcriptome/gencode_transcriptome.fa",
    output:
        gencode_decoys="results/prepare/create_decoys_salmon_index_novel/gencode_decoys.txt",
        gencode_decoyed_transcriptome="results/prepare/create_decoys_salmon_index_novel/gencode_decoyed_sequence.fa",
    params:
        outdir="results/prepare/create_decoys_salmon_index_novel",
    log:
        "logs/prepare/create_decoys_salmon_index_novel/out.log",
    shell:
        """
        mkdir -p {params.outdir};
        grep "^>" {input.gencode_genome} | cut -d " " -f 1 \
            > {output.gencode_decoys} 2>> {log};
        sed -i.bak -e 's/>//g' {output.gencode_decoys} 2>> {log};
        cat {input.gencode_transcriptome} {input.gencode_genome} \
            > {output.gencode_decoyed_transcriptome} 2>> {log}
        """


rule prepare_create_salmon_index_novel:
    input:
        decoys="results/prepare/create_decoys_salmon_index_novel/gencode_decoys.txt",
        decoyed_transcriptome="results/prepare/create_decoys_salmon_index_novel/gencode_decoyed_sequence.fa",
    output:
        directory("results/prepare/create_salmon_index_novel/gencode"),
    params:
        k=config["salmon_index_k"],
    threads: config["align_map_bam_threads"]
    log:
        "logs/prepare/create_salmon_index_novel/gencode.log",
    conda:
        "../envs/standalone/salmon.yaml"
    shell:
        """
        salmon --version > {log};
        salmon index -t {input.decoyed_transcriptome} \
            -i {output} -d {input.decoys} \
            --keepDuplicates \
            -k {params.k} -p {threads}  &>> {log}
        """


rule prepare_create_salmon_index:
    input:
        decoys="results/prepare/create_decoys_salmon_index/{type}_decoys.txt",
        decoyed_transcriptome="results/prepare/create_decoys_salmon_index/{type}_decoyed_sequence.fa",
    output:
        directory("results/prepare/create_salmon_index/{type}"),
    params:
        k=config["salmon_index_k"],
    threads: config["align_map_bam_threads"]
    log:
        "logs/prepare/create_salmon_index/{type}.log",
    conda:
        "../envs/standalone/salmon.yaml"
    shell:
        """
        salmon --version > {log};
        salmon index -t {input.decoyed_transcriptome} \
            -i {output} -d {input.decoys} \
            --keepDuplicates \
            -k {params.k} -p {threads}  &>> {log}
        """


rule prepare_create_refgene_liqa_sirv:
    input:
        "results/prepare/adjust_sirv_names/sirv_set_four.gtf",
    output:
        "results/prepare/create_refgene_liqa_sirv/sirv.refgene",
    params:
        k=config["salmon_index_k"],
    threads: config["align_map_bam_threads"]
    log:
        "logs/prepare/create_refgene_liqa_sirv/out.log",
    benchmark:
        repeat("benchmarks/prepare/create_refgene_liqa_sirv/out.txt", 5)
    conda:
        "../envs/standalone/liqa.yaml"
    shell:
        """
        liqa --version > {log};
        liqa -task refgene -ref {input} -format gtf -out {output} &>> {log}
        """


rule prepare_compile_lr_kallisto:
    output:
        directory("results/prepare/compile_lr_kallisto"),
    log:
        "logs/prepare/compile_lr_kallisto/out.log",
    conda:
        "../envs/standalone/kallisto.yaml"
    shell:
        """
        mkdir {output}; 
        cd {output};
        git clone https://github.com/pachterlab/kallisto &> ../../../{log};
        cd kallisto;
        mkdir build;
        cd build;
        cmake .. -DUSE_HDF5=ON -DMAX_KMER_SIZE=64 &>> ../../../../../{log};
        make &>> ../../../../../{log}
        """


rule prepare_create_lr_kallisto_index:
    input:
        kallisto_binary="results/prepare/compile_lr_kallisto",
        transcriptome="results/prepare/extract_transcriptomes/{type}_transcriptome.fa",
    output:
        f"results/prepare/create_lr_kallisto_index/{{type}}_k-{config['lr_kallisto_index_k']}.idx",
    params:
        k=config["lr_kallisto_index_k"],
    threads: config["align_map_bam_threads"]
    log:
        "logs/prepare/create_lr_kallisto_index/{type}.log",
    benchmark:
        repeat("benchmarks/prepare/create_lr_kallisto_index/{type}.txt", 1)
    shell:
        """
        {input.kallisto_binary}/kallisto/build/src/kallisto index \
            -k {params.k} -t {threads} -i {output} \
            {input.transcriptome} &> {log}
        """


rule prepare_create_lr_kallisto_index_novel:
    input:
        kallisto_binary="results/prepare/compile_lr_kallisto",
        transcriptome="results/prepare/extract_novel_transcriptome/gencode_transcriptome.fa",
    output:
        f"results/prepare/create_lr_kallisto_index_novel/gencode_novel_k-{config['lr_kallisto_index_k']}.idx",
    params:
        k=config["lr_kallisto_index_k"],
    threads: config["align_map_bam_threads"]
    log:
        "logs/prepare/create_lr_kallisto_index_novel/out.log",
    shell:
        """
        {input.kallisto_binary}/kallisto/build/src/kallisto index \
            -k {params.k} -t {threads} -i {output} \
            {input.transcriptome} &> {log}
        """


rule prepare_install_bambu:
    output:
        "results/prepare/install_bambu/done.txt",
    params:
        version=config["bambu_version"],
    log:
        "logs/prepare/install_bambu/out.log",
    conda:
        "../envs/r/bambu.yaml"
    script:
        "../scripts/r/prepare/install_bambu.R"
