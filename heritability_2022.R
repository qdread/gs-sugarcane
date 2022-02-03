# SUGARCANE HERITABILITY
# Author: Quentin Read
# Date: 28 Jan 2022

# Copy the code from GS_sugarcane_2022.R to load phenotype and genotype data
# Then calculate heritability
# Note: clone checks are removed before estimating the variance components, and missing values are excluded (not imputed)

# Load packages -----------------------------------------------------------

library(readxl)
library(data.table)
library(purrr)
library(rrBLUP) 
library(furrr)

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

# One-hot encoding of crop cycle column to use in mixed.solve() as fixed effect design matrix
X <- dcast(data.table(id = 1:nrow(phenotypes), crop_cycle = phenotypes$crop_cycle), id ~ crop_cycle, length)
X[, id := NULL]

# One-hot encoding of clone column to use in mixed.solve()
Z <- dcast(data.table(id = 1:nrow(phenotypes), Clone = phenotypes$Clone), id ~ Clone, length)
Z[, id := NULL]

# Sort columns of Z to match additive matrix calculated later
Z <- Z[, ..clone_columns]

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

# Loop through to calculate h2 --------------------------------------------

# We need to calculate both broad sense and narrow sense heritability for each crop cycle separately and combined

physical_traits <- c("stkwt_kg", "diam", "Brix", "Fiber", "Pol", "Sucrose", "Purity", "stalk_ha")
economic_traits <- c("TCH", "TRS", "CRS", "TSH", "EI")
crop_cycles <- c("Plant Cane_2017", "Ratoon 1_2018", "Ratoon 2_2019")

# Function to get heritability from a model fit
h2 <- function(fit) fit[['Vu']]/(fit[['Vu']] + fit[['Ve']])

# Function to get narrow and broad heritability for a given trait, both with all crop cycles combined and separately
get_all_h2 <- function(trait_to_use) {
  message('Trait: ', trait_to_use)
  
  y <- phenotypes[[trait_to_use]]
  
  # All crop cycles combined in one model. Crop cycle fixed effect, genotype random effect
  # Pass genotype one-hot encoded as random effects Z, crop-cycle one-hot encoded as fixed effects X
  # Pass additive matrix A as random effects covariance matrix K for narrow sense
  # Do not pass anything as K matrix for broad sense. It will be assumed to be the identity matrix.
  fit_combined_narrow <- mixed.solve(y = y, Z = Z, K = A, X = X, method = 'REML')
  fit_combined_broad <- mixed.solve(y = y, Z = Z, X = X, method = 'REML')
  h2_combined <- data.frame(crop_cycle = 'combined', h2_narrow = h2(fit_combined_narrow), H2_broad = h2(fit_combined_broad))
  
  # Each crop cycle as a separate model.
  # Pass genotype as random effects Z and covariance K as above, but do not pass any fixed effects X
  h2_by_cycle <- map_dfr(crop_cycles, function(cycle) {
    idx <- phenotypes$crop_cycle == cycle
    fit_cycle_narrow <- mixed.solve(y = y[idx], Z = Z[idx, ], K = A, method = 'REML')
    fit_cycle_broad <- mixed.solve(y = y[idx], Z = Z[idx, ], method = 'REML')
    data.frame(crop_cycle = cycle, h2_narrow = h2(fit_cycle_narrow), H2_broad = h2(fit_cycle_broad))
  })
  
  data.table(trait = trait_to_use, rbind(h2_combined, h2_by_cycle))
}

# Calculate h2 for all traits in parallel across 4 cores
plan(multisession(workers = 4))

h2_all_traits <- future_map_dfr(c(physical_traits, economic_traits), get_all_h2)
fwrite(h2_all_traits, 'project/output/h2_all_traits.csv')
