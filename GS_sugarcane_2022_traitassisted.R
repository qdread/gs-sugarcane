# SUGARCANE TRAIT-ASSISTED GS
# Author: Quentin Read
# Date: 24 Feb 2022

# The goal of this script is to do trait-assisted genomic selection,
# using a combination of genotype marker data plus traits from previous crop cycle(s)
# to predict traits in a later crop cycle.

# See Fernandes et al. 2018, Theoretical and Applied Genetics

# Load packages -----------------------------------------------------------

library(readxl)
library(data.table)
library(purrr)
library(rrBLUP) 
library(sommer)
library(furrr)
library(lme4)
library(glue)

source('GS_sugarcane_2022_fns.R') # Load needed functions

# Read genotype and phenotype data ----------------------------------------

genotypes <- fread('project/data/sugarcane.10501.SNPs.432.Inds.CloneNames.hmp.txt', sep = '\t')

# Read each sheet from phenotype XLSX and combine to single data frame
# Note there are multiple different NA flags in the data
pheno_file <- 'project/data/Phenotype_updated_2017-18_IL_analysis_120921.xlsx'
sheet_names <- excel_sheets(pheno_file)
phenotype_sheets <- map(sheet_names, ~ cbind(crop_cycle = ., read_xlsx(pheno_file, sheet = ., na = c('.', '-', '-!'))))

# Correct column name, changing Diam to diam (bugfix)
for (i in 1:length(phenotype_sheets)) names(phenotype_sheets[[i]])[names(phenotype_sheets[[i]]) == 'Diam'] <- 'diam'

phenotypes <- rbindlist(phenotype_sheets, fill = TRUE)

# remove clone checks and Ratoonability
check_IDs <- c('CP96-1252', 'CP00-1101')
phenotypes <- phenotypes[!Clone %in% check_IDs & !crop_cycle %in% 'Ratoonability']

# Find columns representing clones
id_columns <- c("rs.", "alleles", "chrom", "pos", "strand", "assembly.", "center", "protLSID", "assayLSID", "panelLSID", "QCcode")
clone_columns <- setdiff(names(genotypes), id_columns)

# Convert SNPs to numerical format ----------------------------------------

# Exclude rows with >2 alleles
genotypes <- genotypes[nchar(alleles) <= 3]
geno_mat_raw <- as.matrix(genotypes[, ..clone_columns])

# Get ID of first and second allele. Map them to 0 and 2 respectively. Map all other values (heterozygotes) within the row to 1.
alleles_list <- strsplit(genotypes$alleles, '/')
allele1 <- map_chr(alleles_list, 1)
allele2 <- map_chr(alleles_list, 2)

geno_mat <- map2(asplit(geno_mat_raw, 1), alleles_list, function(x, allele) fcase(x == allele[1], 0, x == allele[2], 2, default = 1))
geno_mat <- do.call(cbind, geno_mat) # Transpose the matrix so rows are individuals and columns SNPs

# Set row names of genotype matrix to names of clones
dimnames(geno_mat) <- list(clone_columns, NULL)

# Get additive genetic relationship matrix. Recode to -1, 0, 1 format for rrBLUP by subtracting 1.
A <- A.mat(geno_mat - 1, n.core=3, shrink=FALSE, return.imputed=FALSE)

# Get BLUPs for each trait and crop cycle ---------------------------------

# Group columns by ID, physical traits, and economic traits
physical_traits <- c("stkwt_kg", "diam", "Brix", "Fiber", "Pol", "Sucrose", "Purity", "stalk_ha")
economic_traits <- c("TCH", "TRS", "CRS", "TSH", "EI")

phenotypes_long <- phenotypes[!crop_cycle %in% 'Ratoonability', mget(c("crop_cycle", "Clone", "Rep", "Row", "Column", physical_traits, economic_traits))] |>
  melt(id.vars = c("crop_cycle", "Clone", "Rep", "Row", "Column"), variable.name = 'trait')

# Get the BLUPs!
pheno_blups <- phenotypes_long[, blup_trait(.SD, crop_cycle_fixef = FALSE), by = .(trait, crop_cycle)]
setnames(pheno_blups, c('trait', 'crop_cycle', 'Clone', 'BLUP'))

# Remove the clone checks 
check_IDs <- c('CP96-1252', 'CP00-1101')
pheno_blups <- pheno_blups[!Clone %in% check_IDs]

# Change crop cycle names to the short name
pheno_blups[, crop_cycle := fcase(
  crop_cycle == 'Plant Cane_2017', 'PlantCane',
  crop_cycle == 'Ratoon 1_2018', 'Ratoon1',
  crop_cycle == 'Ratoon 2_2019', 'Ratoon2'
)]

# Do the trait-assisted GS ------------------------------------------------

n_iter <- 25 
n_folds <- 5

# Repeat for each:
# - crop cycle (predict 1st ratoon with plant cane, predict 2nd ratoon with combination of plant cane and 1st ratoon)
# - trait
# - model

# Then in each case do n_iter iterations, and in each iteration the n_folds CV folds.

combos <- CJ(iter = 1:n_iter, trait = c(physical_traits, economic_traits), crop_cycle = c('Ratoon1', 'Ratoon2'))

# Do the GS. Write observed and predicted phenotypes and prediction accuracy metrics with each iteration.
# If the prediction metric file already exists, skip that iteration (this allows the script to be rerun)
traitgs_pred_metrics <- future_pmap(combos, function(iter, trait, crop_cycle) {
  metric_file_name <- glue('project/output/TAGS/metrics_{trait}_{crop_cycle}_{iter}.csv')
  if (!file.exists(metric_file_name)) {
    pred_vals <- trait_assisted_gs(
      GD = geno_mat, PD = pheno_blups, crop_cycle_to_use = crop_cycle, trait_to_use = trait, k = n_folds
    )
    fwrite(pred_vals, glue('project/output/TAGS/phenotypes_{trait}_{crop_cycle}_{iter}.csv'))
    pred_metrics <- with(pred_vals, calc_metrics(Y_obs, Y_pred))
    fwrite(pred_metrics, metric_file_name)
  } else {
    pred_metrics <- fread(metric_file_name)
  }
  return(pred_metrics)
}, .options = furrr_options(seed = 313))

combos[, metrics := traitgs_pred_metrics]
results <- unnest_dt(combos, col = metrics, id = .(iter, trait, crop_cycle))

fwrite(results, 'project/output/TAGS_all_metrics.csv')

