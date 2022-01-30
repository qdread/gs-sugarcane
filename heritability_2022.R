# SUGARCANE HERITABILITY
# Author: Quentin Read
# Date: 28 Jan 2022

# Copy the code from GS_sugarcane_2022.R to load phenotype and genotype data
# Then calculate heritability

# Load packages -----------------------------------------------------------

library(readxl)
library(data.table)
library(purrr)
library(sommer) 

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


# Find broad-sense h2 -----------------------------------------------------


# Find narrow-sense h2 ----------------------------------------------------




