configfile: "config/config.yaml"


singularity: f"docker://condaforge/mambaforge:{config['mambaforge_version']}"


rule annotate_run_sqanti:
    input:
        query="results/discover/fix_bambu_gene_ids/transcriptome.gtf",
        reference="results/prepare/adjust_transcriptome_assembly_names/gencode.v45.primary_assembly.annotation.named.gtf",
        genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        polya_motifs="results/prepare/download_sqanti_polya_motif_list/motifs.txt",
    output:
        classification="results/annotate/run_sqanti/kinnex_wtc_11_classification.txt",
        junctions="results/annotate/run_sqanti/kinnex_wtc_11_junctions.txt",
    params:
        output_folder="results/annotate/run_sqanti/kinnex_wtc_11",
        sample_id="kinnex_wtc_11",
    threads: config["parallel_threads"]
    log:
        "logs/annotate/run_sqanti/out.log",
    conda:
        "../envs/py/sqanti.yaml"
    shell:
        """
        conda list &> {log};
        python workflow/scripts/py/sqanti3_wrapper.py \
                {input.query} {input.reference} {input.genome} \
                --force_id_ignore \
                --skipORF \
                --polyA_motif_list {input.polya_motifs} \
                -o {params.sample_id} -d {params.output_folder} \
                --cpus {threads} --report both \
                --min_ref_len 0 &>> {log}
                """


rule annotate_run_orfanage:
    input:
        query_transcriptome="results/discover/fix_bambu_gene_ids/transcriptome.gtf",
        genome="results/prepare/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        reference_transcriptome="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
    output:
        f"results/annotate/run_orfanage/transcriptome.gtf",
    params:
        mode=config["call_orfs_orfanage_mode"],
    threads: config["parallel_threads"]
    log:
        "logs/call_orfs/run_orfanage/transcriptome.log",
    conda:
        "../envs/standalone/orfanage.yaml"
    shell:
        """
        conda list &> {log};
        orfanage \
            --mode {params.mode} \
            --reference {input.genome} \
            --query {input.query_transcriptome} \
            --output {output} \
            --threads {threads} \
            --keep_cds \
            {input.reference_transcriptome} \
            &>> {log}
        """


rule annotate_prepare_sqanti_protein_gffread:
    input:
        sample="results/annotate/run_orfanage/transcriptome.gtf",
        reference="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
    output:
        sample="results/annotate/prepare_sqanti_protein_gffread/transcriptome.gtf",
        reference="results/annotate/prepare_sqanti_protein_gffread/gencode.gtf",
    log:
        "logs/annotate/prepare_sqanti_protein_gffread/out.log",
    conda:
        "../envs/standalone/gffread.yaml"
    shell:
        """
        conda list &> {log};
        gffread -C -E -T {input.sample} -o {output.sample} &>> {log};
        gffread -C -E -T {input.reference} -o {output.reference} &>> {log}
        """


rule annotate_prepare_sqanti_protein:
    input:
        query="results/annotate/prepare_sqanti_protein_gffread/transcriptome.gtf",
        reference="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
    output:
        query="results/annotate/prepare_sqanti_protein/transcriptome.gtf",
        reference="results/annotate/prepare_sqanti_protein/gencode.gtf",
    log:
        "logs/annotate/prepare_sqanti_protein/out.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/prepare_sqanti_protein.R"


rule annotate_prepare_sqanti_protein_construct_gene_features:
    input:
        query="results/annotate/prepare_sqanti_protein/transcriptome.gtf",
        reference="results/annotate/prepare_sqanti_protein/gencode.gtf",
    output:
        query="results/annotate/prepare_sqanti_protein_construct_gene_features/transcriptome.gtf",
        reference="results/annotate/prepare_sqanti_protein_construct_gene_features/gencode.gtf",
    log:
        "logs/annotate/prepare_sqanti_protein_construct_gene_features/out.log",
    conda:
        "../envs/standalone/gffread.yaml"
    shell:
        """
        conda list &> {log};
        gffread -E --keep-genes {input.query} -o- 2>> {log} | \
            gffread -E - -T -o {output.query} &>> {log};
        gffread -E --keep-genes {input.reference} -o- 2>> {log} | \
            gffread -E - -T -o {output.reference} &>> {log}
        """


# Adapted from: https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/460d53ae716419438b31a3a8052f3fce2e6175ce/main.nf
rule annotate_run_sqanti_protein:
    input:
        query="results/annotate/prepare_sqanti_protein_gffread/transcriptome.gtf",
        reference="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        query_cds="results/annotate/prepare_sqanti_protein_construct_gene_features/gencode.gtf",
        reference_cds="results/annotate/prepare_sqanti_protein_construct_gene_features/transcriptome.gtf",
    output:
        "results/annotate/run_sqanti_protein/kinnex_wtc_11.sqanti_protein_classification.tsv",
        outdir=directory("results/annotate/run_sqanti_protein"),
    params:
        sample_id="kinnex_wtc_11",
    log:
        f"logs/annotate/run_sqanti_protein/out.log",
    conda:
        "../envs/py/sqanti.yaml"
    shell:
        """
        conda list &> {log};
        python workflow/scripts/py/sqanti3_protein.py \
                {input.query} \
                {input.query_cds} \
                {input.reference} \
                {input.reference_cds} \
                -d {output.outdir} \
                -p {params.sample_id} &>> {log}
        """


rule annotate_postprocess_sqanti_protein:
    input:
        protein_classification="results/annotate/run_sqanti_protein/kinnex_wtc_11.sqanti_protein_classification.tsv",
        query="results/annotate/prepare_sqanti_protein_gffread/transcriptome.gtf",
        reference="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
    output:
        "results/annotate/postprocess_sqanti_protein/transcriptome.protein_classification.tsv",
        outdir=directory("results/annotate/postprocess_sqanti_protein/"),
    params:
        sample_id="kinnex_wtc_11",
    log:
        "logs/annotate/postprocess_sqanti_protein/out.log",
    conda:
        "../envs/py/sqanti_postprocess.yaml"
    shell:
        """
        conda list &> {log};
        python workflow/scripts/py/get_gc_exon_and_5utr_info.py \
                --gencode_gtf {input.reference} \
                --odir {output.outdir} &>> {log};
        python workflow/scripts/py/classify_5utr_status.py \
                --gencode_exons_bed {output.outdir}/gencode_exons_for_cds_containing_ensts.bed \
                --gencode_exons_chain {output.outdir}/gc_exon_chain_strings_for_cds_containing_transcripts.tsv \
                --sample_cds_gtf {input.query} \
                --odir {output.outdir} &>> {log};
        python workflow/scripts/py/merge_5utr_info_to_pclass_table.py \
            --name {params.sample_id} \
            --utr_info {output.outdir}/pb_5utr_categories.tsv \
            --sqanti_protein_classification {input.protein_classification} \
            --odir {output.outdir} &>> {log};
        python workflow/scripts/py/protein_classification.py \
            --sqanti_protein {output.outdir}/{params.sample_id}.sqanti_protein_classification_w_5utr_info.tsv \
            --name {params.sample_id} \
            --dest_dir {output.outdir} &>> {log}
        """
