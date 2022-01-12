
# Load packages -----------------------------------------------------------

library(readxl)
library(data.table)
library(purrr)
library(rrBLUP)
library(kernlab)
library(e1071)
library(randomForest)
library(BGLR)
library(sommer) # Note: some functions in rrBLUP have the same name as some in this package.

source('GS_sugarcane_2022_fns.R')

# Read genotype and phenotype data ----------------------------------------

genotypes <- fread('project/data/sugarcane.10501.SNPs.432.Inds.CloneNames.hmp.txt', sep = '\t')

# Read each sheet from phenotype XLSX and combine to single data frame
# Note there are multiple different NA flags in the data
sheet_names <- excel_sheets('project/data/Phenotype_updated_2017-18_IL_analysis_120921.xlsx')
phenotype_sheets <- map(sheet_names, ~ cbind(crop_cycle = ., read_xlsx('project/data/Phenotype_updated_2017-18_IL_analysis_120921.xlsx', sheet = ., na = c('.', '-', '-!'))))
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

# Average traits by clone -------------------------------------------------

# Group columns by ID, physical traits, and economic traits
phen_id_columns <- c("variable", "Clone", "Rep", "Tier", "Plot", "Crop", "Row", "Column")
physical_traits <- c("stkwt_kg", "diam", "Brix", "Fiber", "Pol", "Sucrose", "Purity", "stalk_ha")
economic_traits <- c("TCH", "TRS", "CRS", "TSH", "EI")

pheno_means <- phenotypes[, lapply(.SD, mean, na.rm = TRUE), by = .(crop_cycle, Clone), .SDcols = c(physical_traits, economic_traits)]

# Remove the clone checks 
check_IDs <- c('CP96-1252', 'CP00-1101')
pheno_means <- pheno_means[!is.na(Clone) & !Clone %in% check_IDs]

# Setup data structures for output ----------------------------------------

n_iter <- 25
n_folds <- 5

# Repeat for each:
# - crop cycle (plant cane, 1st ratoon, 2nd ratoon)
# - trait
# - training set size c(0.20, 0.30, 0.50, 0.60, 0.80, 0.90). Instead for now I will just do 5fold CV? 
# - marker density c(0.20, 0.30, 0.50, 0.60, 0.80, 0.90, 1)
# - models (5 models, same as those done by Marcus)

# Then in each case do 25 iterations.

testgs <- gs_all(GD = geno_mat, PD = pheno_means, 
                 crop_cycle_to_use = 'Plant Cane_2017', trait = 'stkwt_kg', k = n_folds, rand_seed = 333, marker_density = 1)
