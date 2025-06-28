# Quantitative isoform profiling using deep coverage long-read RNA sequencing across early endothelial differentiation

## Abstract

PacBio long-read RNA sequencing resolves transcripts with greater clarity than short-read technologies, yet its quantitative performance remains under-evaluated at scale. Here, we benchmark the high-throughput PacBio Kinnex platform against Illumina short-read RNA-seq using matched, deeply sequenced datasets across a time course of endothelial cell differentiation. Compared to Illumina, Kinnex achieved comparable gene-level quantification and more accurate transcript discovery and transcript quantification. While Illumina detected more transcripts overall, many reflected potentially unstable or ambiguous estimates in complex genes. Kinnex largely avoids these issues, producing more reliable differential transcript expression (DTE) calls, despite a mild bias against short transcripts (<1.25 kb). When correcting Illumina for inferential variability, Kinnex and Illumina quantifications were highly concordant, demonstrating equivalent performance. We also benchmarked long-read tools, nominating Oarfish as the most efficient for our Kinnex data. Together, our results establish Kinnex as a reliable platform for full-length transcript quantification

## Reproducibility

To reproduce our results, you need to have `snakemake>=7.30.1`, `singularity-ce` or `apptainer` and a recent `conda` (or alternative frontend, such as `mamba`) version installed. Starting from the raw data (which is controlled access, see `Data availibility`) you can easily reproduce all of our analyses, including the tables and figures contained in the paper by running Snakemake:

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

### Code

In addition to here on Github, a copy of this repo is also available from [Zenodo](https://doi.org/10.5281/zenodo.15502054/).

## Citation

Our manuscript is still under review. For now, you may find a BibTex entry for our preprint below.

```tex
@article{wissel2025systematic,
  title={A Systematic Benchmark of High-Accuracy PacBio Long-Read RNA Sequencing for Transcript-Level Quantification},
  author={Wissel, David and Mehlferber, Madison M and Nguyen, Khue M and Pavelko, Vasilii and Tseng, Elizabeth and Robinson, Mark D and Sheynkman, Gloria M},
  journal={bioRxiv},
  pages={2025--06},
  year={2025},
  publisher={Cold Spring Harbor Laboratory}
}
```

## Contact

In case of any questions, please reach out to [Madison](mm5db@virginia.edu) and [David](dwissel@ethz.ch) or open an issue in this repo.
