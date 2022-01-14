# GENOMIC SELECTION SUGARCANE: SCRIPT TO RUN GS MODELS
# Author: Quentin Read
# Date: 12 January 2022
# ----------------------------------------------------

# Load packages -----------------------------------------------------------

library(readxl)
library(data.table)
library(purrr)
library(furrr)
library(glue)
library(lme4)
library(uuid)
library(rrBLUP)
library(kernlab)
library(e1071)
library(randomForest)
library(BGLR)
library(sommer) # Note: some functions in rrBLUP have the same name as some in this package.

source('GS_sugarcane_2022_fns.R') # Load needed functions

# Set up parallel processing
n_tasks <- as.numeric(Sys.getenv('SLURM_NTASKS'))
options(mc.cores = n_tasks)
plan(multicore(workers = n_tasks))

# Read genotype and phenotype data ----------------------------------------

genotypes <- fread('project/data/sugarcane.10501.SNPs.432.Inds.CloneNames.hmp.txt', sep = '\t')

# Read each sheet from phenotype XLSX and combine to single data frame
# Note there are multiple different NA flags in the data
pheno_file <- 'project/data/Phenotype_updated_2017-18_IL_analysis_120921.xlsx'
sheet_names <- excel_sheets(pheno_file)
phenotype_sheets <- map(sheet_names, ~ cbind(crop_cycle = ., read_xlsx(pheno_file, sheet = ., na = c('.', '-', '-!'))))
phenotypes <- rbindlist(phenotype_sheets, fill = TRUE)

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


# Get BLUPs for each trait and crop cycle ---------------------------------

# Group columns by ID, physical traits, and economic traits
physical_traits <- c("stkwt_kg", "diam", "Brix", "Fiber", "Pol", "Sucrose", "Purity", "stalk_ha")
economic_traits <- c("TCH", "TRS", "CRS", "TSH", "EI")

phenotypes_long <- phenotypes[!crop_cycle %in% 'Ratoonability', mget(c("crop_cycle", "Clone", "Rep", "Row", "Column", physical_traits, economic_traits))] |>
  melt(id.vars = c("crop_cycle", "Clone", "Rep", "Row", "Column"), variable.name = 'trait')

# Function to calculate BLUPs for each trait. If only one crop cycle has that trait, do not fit crop cycle in the model.
blup_trait <- function(dat) {
  if (length(unique(dat$crop_cycle)) > 1) {
    lmm <- lmer(value ~ 0 + crop_cycle + (1|Clone) + (1|Row) + (1|Column), data = dat,
                control = lmerControl(optimizer = 'bobyqa'))
    blup <- outer(ranef(lmm)[['Clone']][['(Intercept)']], fixef(lmm), `+`)
    data.frame(Clone = row.names(ranef(lmm)[['Clone']]), blup)
  } else {
    lmm <- lmer(value ~ 1 + (1|Clone) + (1|Row) + (1|Column), data = dat,
                control = lmerControl(optimizer = 'bobyqa'))
    blup <- outer(ranef(lmm)[['Clone']][['(Intercept)']], fixef(lmm), `+`)
    data.frame(Clone = row.names(ranef(lmm)[['Clone']]), blup, as.numeric(NA), as.numeric(NA))
  }
}

# First and second ratoons do not have diameter data. Remove those rows
phenotypes_long <- phenotypes_long[!(trait %in% 'diam' & crop_cycle %in% c('Ratoon 1_2018', 'Ratoon 2_2019'))]

# Get the BLUPs!
pheno_blups <- phenotypes_long[, blup_trait(.SD), by = trait]
setnames(pheno_blups, c('trait', 'Clone', 'PlantCane', 'Ratoon1', 'Ratoon2'))

# Remove the clone checks 
check_IDs <- c('CP96-1252', 'CP00-1101')
pheno_blups <- pheno_blups[!Clone %in% check_IDs]


# Run GS models -----------------------------------------------------------

n_iter <- 25
n_folds <- 5

# Repeat for each:
# - crop cycle (plant cane, 1st ratoon, 2nd ratoon)
# - trait
# - training set size c(0.20, 0.30, 0.50, 0.60, 0.80, 0.90). Instead for now I will just do 5fold CV? 
# - marker density c(0.20, 0.30, 0.50, 0.60, 0.80, 0.90, 1)
# - models (5 models, same as those done by Marcus)

# Then in each case do 25 iterations.

# For now, just do iterations, models, crop cycles, and traits. Do not vary training size or marker density.

combos <- CJ(iter = 1:n_iter, trait = c(physical_traits, economic_traits), crop_cycle = c('PlantCane', 'Ratoon1', 'Ratoon2'))

# 1st and 2nd ratoons do not have diameter data so remove those from the combinations
combos <- combos[!(trait %in% 'diam' & crop_cycle %in% c('Ratoon_1', 'Ratoon_2'))]

# Do the GS. Write observed and predicted phenotypes and prediction accuracy metrics with each iteration.
gs_pred_metrics <- future_pmap(combos, function(iter, trait, crop_cycle) {
  pred_vals <- gs_all(GD = geno_mat, PD = pheno_blups, 
                      crop_cycle_to_use = crop_cycle, trait_to_use = trait, k = n_folds, marker_density = 1)
  fwrite(pred_vals, glue('project/output/phenotypes_{trait}_{crop_cycle}_{iter}.csv'))
  pred_metrics <- pred_vals[, calc_metrics(Y_obs, Y_pred), by = model]
  fwrite(pred_metrics, glue('project/output/metrics_{trait}_{crop_cycle}_{iter}.csv'))
  return(pred_metrics)
}, .options = furrr_options(seed = 777))

combos[, metrics := gs_pred_metrics]
results <- unnest_dt(combos, col = metrics, id = .(iter, trait, crop_cycle))

fwrite(results, 'project/output/all_metrics.csv')
