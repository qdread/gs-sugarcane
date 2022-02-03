# GENOMIC SELECTION SUGARCANE: Combine results and write to single file to export
# Author: Quentin Read
# Date: 13 January 2022
# -------------------------------------------------------------------------------

library(data.table)
library(glue)

physical_traits <- c("stkwt_kg", "diam", "Brix", "Fiber", "Pol", "Sucrose", "Purity", "stalk_ha")
economic_traits <- c("TCH", "TRS", "CRS", "TSH", "EI")

combos <- CJ(iter = 1:25, trait = c(physical_traits, economic_traits), crop_cycle = c('PlantCane', 'Ratoon1', 'Ratoon2'))
combos[, file_name := with(combos, glue('project/output/metrics_{trait}_{crop_cycle}_{iter}.csv'))]
combos <- combos[file.exists(file_name)]
metrics <- combos[, fread(file_name), by = .(iter, trait, crop_cycle)]

fwrite(metrics, 'project/metrics_18jan2022.csv')


# Export observed and predicted values ------------------------------------

library(data.table)
library(glue)

physical_traits <- c("stkwt_kg", "diam", "Brix", "Fiber", "Pol", "Sucrose", "Purity", "stalk_ha")
economic_traits <- c("TCH", "TRS", "CRS", "TSH", "EI")

combos <- CJ(iter = 1:25, trait = c(physical_traits, economic_traits), crop_cycle = c('PlantCane', 'Ratoon1', 'Ratoon2'))
combos[, file_name := with(combos, glue('project/output/phenotypes_{trait}_{crop_cycle}_{iter}.csv'))]
combos <- combos[file.exists(file_name)]
phenotypes <- combos[, fread(file_name), by = .(iter, trait, crop_cycle)]
phenotypes <- phenotypes[, .(trait, crop_cycle, iter, model, fold, Clone, Y_obs, Y_pred)]
setnames(phenotypes, old = 'iter', new = 'iteration')

fwrite(phenotypes, 'project/observed_predicted_phenotypes_03feb2022.csv')