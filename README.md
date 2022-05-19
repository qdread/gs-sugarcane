# Accuracy of genomic prediction of yield and sugar traits in sugarcane

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6563705.svg)](https://doi.org/10.5281/zenodo.6563705)

This repository accompanies the manuscript:

Islam, Md., P. H. McCord, Q. D. Read, L. Qin, A. Lipka, S. Sood, J. Todd, and M. Olatoye. Accuracy of genomic prediction of yield and sugar traits in *Saccharum* spp. hybrids. *Agriculture*, in preparation.

## Project contributors

- Md. Sariful Islam, USDA ARS (PI)
- Per McCord, Washington State University
- Quentin D. Read, USDA ARS/NC State University (analyst)
- Lifang Qin, USDA ARS
- Alexander Lipka, University of Illinois
- Sushma Sood, USDA ARS
- James Todd, USDA ARS
- Marcus Olatoye, University of Illinois


## Data

All code and data required to reproduce the analyses in the manuscript are included here in this repository (the input data are only ~10 MB total). 

Data are found in the `project/data` folder in this repo. There are two files:

- `Phenotype_updated_2017-18_IL_analysis_120921.xlsx`: Phenotype data in a MS Excel spreadsheet with columns of trait values from each crop cycle in a separate tab. This also includes experimental layout information (the replicate, plot, row, and column for each measured individual).
- `sugarcane.10501.SNPs.432.Inds.CloneNames.hmp.txt`: Genotype data in a single tab-separated text file with one column for each of the 432 clones and one row for each SNP locus.

## How to run the analysis

### Model fitting

There are five major steps to the analysis that can be run in any order. The first three scripts are set up to be run on a remote Slurm cluster over 4 nodes with 36 processor cores each using the `rslurm` package, the fourth script submits the job to the cluster directly without using `rslurm`, and the last script is run locally assuming 4 processor cores are available.

- `GS_sugarcane_2022_runGS_rslurm.R`: fits all genomic selection models for all combinations of crop cycle and trait and performs 5-fold cross-validation to generate out-of-sample predictions. The entire cross-validation procedure is repeated 25 times for every combination of crop cycle and trait.
- `GS_sugarcane_2022_markerdensity_rslurm.R`: repeats the same cross-validation procedure (k = 5 with 25 iterations) for every combination of crop cycle and trait, but the genotypic marker dataset is subsampled to 20%, 30%, 50%, 60%, 80%, and 90% of its original size.
- `GS_sugarcane_2022_trainingsize_rslurm.R`: repeats the same cross-validation procedure (k = 5 with 25 iterations) for every combination of crop cycle and trait, but the training dataset of genotypes is subsampled to 20%, 30%, 50%, 60%, 80%, and 90% of its original size within each cross-validation fold (i.e., the 80% of the data used for fitting the model within each fold are further subsampled).
- `GS_sugarcane_2022_traitassisted.R`: performs trait-assisted genomic selection for each trait for the two later crop cycles. For the first ratoon, the traits measured at plant cane stage are included in the multivariate response, and for the second ratoon, both the plant cane and first ratoon traits are included. Submit this using the job submission script `doalltags.sh`.
- `heritability_2022.R`: calculates broad-sense and narrow-sense heritabilities for each trait for the three crop cycles separately, and for all crop cycles combined.

Each of the four scripts above sources two other scripts:

- `GS_sugarcane_2022_fns.R`: defines all functions needed to fit the models.
- `GS_sugarcane_2022_loaddata_rslurm.R`: loads all necessary R packages and data, and pre-processes the data.

### Visualizing results

After running the model fitting scripts above, run the following script:

- `export_all_metrics.R`: combines the many temporary output files from each iteration of the analyses above into a single CSV for each main component of the analysis.

Results are visualized and described, using figures, tables, and text, in three RMarkdown notebooks which can be rendered into HTML documents once the results have been generated:

- `results_GS.Rmd`: main genomic selection results
- `results_TS_MD.Rmd`: results of the model fitting varying training size and marker density
- `results_TAGS.Rmd`: trait-assisted genomic selection results

PNG figures for the manuscript are created in the following scripts:

- `figs_for_ms.R`: the majority of the figures, comparing model performance across models and traits
- `figs_for_ms_density_scatter.R`: density plots of traits and scatterplots comparing model performance for each trait with heritability and trait variation

*This document last modified by QDR, 07 April 2022*.
