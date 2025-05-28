# Quantitative isoform profiling using deep coverage long-read RNA sequencing across early endothelial differentiation

## Abstract

Long-read RNA-sequencing enables the profiling of full-length transcripts, but its quantification accuracy data has still not been robustly established. This is especially true for PacBio lrRNA-seq data, which were previously only available at low to moderate depth. Using a high-depth PacBio Kinnex lrRNA- seq dataset, sample-matched with Illumina short-read RNA-seq, we performed rigorous benchmarking to characterize quantification accuracy between platforms on a dataset representing differentiation of induced pluripotent stem cells into primordial endothelial cells. We identified biases impacting transcript quantification, including inferential variability within Illumina data, which can bias transcript abundance estimates for genes with complex splicing, as well as length biases in Kinnex data. Overall, PacBio and Illumina quantifications were strongly concordant, supporting that PacBio Kinnex is a reliable method for transcriptome profiling and enabling downstream biological analyses.

## Reproducibility

To reproduce our results, you need to have `snakemake>=7.30.1` and a recent `singularity-ce` or `apptainer` version installed.

```sh
snakemake --use-conda --use-apptainer --cores 12
```

## Data availibility

### Raw data

Raw data is controlled access. For access to the raw sequencing files for both Illumina and Kinnex, please refer to [iTHRIV](https://portal.ithriv.org/#/public_commons/project/25d25b41-87ca-4e41-9573-1ee9375daf1a/dataset/5bcbeb1e-a0e6-4dd3-91ed-107eba1d3061).

### Intermediate data (e.g., count tables, discovered transcriptome, ...)

Intermediate data are available from [Zenodo](https://doi.org/10.5281/zenodo.15502054/). Due to size constraints, not all intermediate data are on Zenodo - if something is missing that you might need, please do not hesitate to contact us.

### Figures and tables

Figures and tables are available from [Zenodo](https://doi.org/10.5281/zenodo.15502054/).

## Citation

Our manuscript is still under review.

## Contact

In case of any questions, please reach out to [Madison](mm5db@virginia.edu) and [David](dwissel@ethz.ch) or open an issue in this repo.
