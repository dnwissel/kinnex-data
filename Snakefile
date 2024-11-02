from snakemake.utils import min_version


configfile: "config/config.yaml"


min_version(config["snakemake_min_version"])


singularity: f"docker://condaforge/mambaforge:{config['mambaforge_version']}"


include: "workflow/rules/prepare.smk"
include: "workflow/rules/align.smk"
include: "workflow/rules/discover.smk"
include: "workflow/rules/downsample.smk"
include: "workflow/rules/quantify_subsampled.smk"


rule all:
    input:
        f"results/prepare/create_lr_kallisto_index/sirv_k-{config['lr_kallisto_index_k']}.idx",
        f"results/prepare/create_lr_kallisto_index/gencode_k-{config['lr_kallisto_index_k']}.idx",
        "results/prepare/create_refgene_liqa_sirv/sirv.refgene",
        "results/prepare/create_salmon_index/sirv",
        "results/prepare/create_salmon_index/gencode",
        "results/prepare/index_star",
        expand(
            "results/prepare/run_trim_galore/{sample}-r1_val_1.fq.gz",
            sample=config["sample_names"],
        ),
        expand(
            "results/prepare/convert_mapped_bams_to_fastq/{sample}/gencode.reads.fastq.gz",
            sample=config["sample_names"],
        ),
        expand(
            "results/prepare/convert_mapped_bams_to_fastq/{sample}/sirv.reads.fastq.gz",
            sample=config["sample_names"],
        ),
        expand(
            "results/align/run_minimap2_{data_type}/{sample}/{sample}.aligned.sorted.bam",
            data_type=["gencode", "sirv"],
            sample=config["sample_names"],
        ),
        expand(
            "results/align/run_minimap2_transcriptome_{data_type}/{sample}/{sample}.aligned.bam",
            data_type=["gencode", "sirv"],
            sample=config["sample_names"],
        ),
        expand(
            "results/align/convert_mapped_bams_to_fastq/{sample}/{data_type}/{sample}.reads_1.fastq.gz",
            data_type=["gencode", "sirv"],
            sample=config["sample_names"],
        ),
        "results/discover/run_bambu/transcriptome.gtf",
        "results/discover/merge_stringtie2/transcriptome.gtf",
        expand(
            "results/quantify_subsampled/run_{method}_subsampled/gencode/{subsample_number}/{read_number}/{sample}/{sample}.quant",
            sample=config["sample_names"],
            method=["oarfish"],
            subsample_number=config["subsample_numbers_gencode"],
            read_number=config["subsample_read_numbers_gencode"],
        ),
        expand(
            "results/quantify_subsampled/run_{method}_subsampled/sirv/{subsample_number}/{read_number}/{sample}/{sample}.quant",
            sample=config["sample_names"],
            method=["oarfish"],
            subsample_number=config["subsample_numbers_sirv"],
            read_number=config["subsample_read_numbers_sirv"],
        ),
        expand(
            "results/quantify_subsampled/run_{method}_subsampled/gencode/{subsample_number}/{read_number}/counts_transcript.txt",
            sample=config["sample_names"],
            method=["bambu"],
            subsample_number=config["subsample_numbers_gencode"],
            read_number=config["subsample_read_numbers_gencode"],
        ),
        expand(
            "results/quantify_subsampled/run_{method}_subsampled/sirv/{subsample_number}/{read_number}/counts_transcript.txt",
            method=["bambu"],
            subsample_number=config["subsample_numbers_sirv"],
            read_number=config["subsample_read_numbers_sirv"],
        ),
        expand(
            "results/quantify_subsampled/run_liqa_subsampled/sirv/{subsample_number}/{read_number}/{sample}.tsv",
            sample=config["sample_names"],
            method=["liqa"],
            subsample_number=config["subsample_numbers_sirv"],
            read_number=config["subsample_read_numbers_sirv"],
        ),
        expand(
            "results/quantify_subsampled/run_{method}_subsampled/sirv/{subsample_number}/{read_number}",
            sample=config["sample_names"],
            method=["isoquant"],
            subsample_number=config["subsample_numbers_sirv"],
            read_number=config["subsample_read_numbers_sirv"],
        ),
        expand(
            "results/quantify_subsampled/run_{method}_subsampled/gencode/{subsample_number}/{read_number}",
            sample=config["sample_names"],
            method=["isoquant"],
            subsample_number=config["subsample_numbers_gencode"],
            read_number=config["subsample_read_numbers_gencode"],
        ),
        expand(
            "results/quantify_subsampled/run_kallisto_long/sirv/{subsample_number}/{read_number}/{sample}",
            sample=config["sample_names"],
            method=["kallisto_long"],
            subsample_number=config["subsample_numbers_sirv"],
            read_number=config["subsample_read_numbers_sirv"],
        ),
        expand(
            "results/quantify_subsampled/run_kallisto_long/gencode/{subsample_number}/{read_number}/{sample}",
            sample=config["sample_names"],
            method=["kallisto_long"],
            subsample_number=config["subsample_numbers_gencode"],
            read_number=config["subsample_read_numbers_gencode"],
        ),
        expand(
            "results/quantify_subsampled/run_{method}_subsampled/gencode/{subsample_number}/{read_number}/{sample}/quant.sf",
            sample=config["sample_names"],
            method=["salmon_illumina", "salmon"],
            subsample_number=config["subsample_numbers_gencode"],
            read_number=config["subsample_read_numbers_gencode"],
        ),
        expand(
            "results/quantify_subsampled/run_{method}_subsampled/sirv/{subsample_number}/{read_number}/{sample}/quant.sf",
            sample=config["sample_names"],
            method=["salmon_illumina", "salmon"],
            subsample_number=config["subsample_numbers_sirv"],
            read_number=config["subsample_read_numbers_sirv"],
        ),
        gencode_transcriptome="results/prepare/extract_transcriptomes/gencode_transcriptome.fa",
        joint_transcriptome="results/prepare/extract_transcriptomes/joint_transcriptome.fa",
        sirv_transcriptome="results/prepare/extract_transcriptomes/sirv_transcriptome.fa",
        gencode_bed="results/prepare/convert_gtfs_to_beds/gencode.bed",
        sirv_bed="results/prepare/convert_gtfs_to_beds/sirv.bed",
