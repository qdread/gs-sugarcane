# GENOMIC SELECTION SUGARCANE: LOAD AND PROCESS DATA
# QDR / 11 Feb 2022
# In this script we load all packages and data, source functions, 
# convert raw marker data to 0,1,2 format, and get BLUPs for phenotypes

# Load packages -----------------------------------------------------------

library(readxl)
library(data.table)
library(purrr)
library(furrr)
library(glue)
library(lme4)
library(Matrix)
library(uuid)
library(rrBLUP)
library(kernlab)
library(e1071)
library(randomForest)
library(BGLR)
library(sommer) # Note: some functions in rrBLUP have the same name as some in this package.

source('GS_sugarcane_2022_fns.R') # Load needed functions

# Set up parallel processing if this code is running on ceres
if (grepl('ceres', system2('hostname', stdout = TRUE))) {
  plan(multicore(workers = length(availableWorkers())))
}

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

# Get the BLUPs!
pheno_blups <- phenotypes_long[, blup_trait(.SD), by = trait]
setnames(pheno_blups, c('trait', 'Clone', 'PlantCane', 'Ratoon1', 'Ratoon2'))

# Remove the clone checks 
check_IDs <- c('CP96-1252', 'CP00-1101')
pheno_blups <- pheno_blups[!Clone %in% check_IDs]