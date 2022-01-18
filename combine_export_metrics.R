# GENOMIC SELECTION SUGARCANE: Combine results and write to single file to export
# Author: Quentin Read
# Date: 13 January 2022
# -------------------------------------------------------------------------------

library(data.table)
library(glue)

source('GS_sugarcane_2022_fns.R') 

physical_traits <- c("stkwt_kg", "diam", "Brix", "Fiber", "Pol", "Sucrose", "Purity", "stalk_ha")
economic_traits <- c("TCH", "TRS", "CRS", "TSH", "EI")

metrics <- CJ(iter = 1:25, trait = c(physical_traits, economic_traits), crop_cycle = c('PlantCane', 'Ratoon1', 'Ratoon2'))
metrics[, file_name := with(metrics, glue('project/output/metrics_{trait}_{crop_cycle}_{iter}.csv'))]
metrics <- metrics[file.exists(file_name)]
metrics[, data := lapply(file_name, fread)]
metrics <- unnest_dt(metrics, col = data, id = .(iter, trait, crop_cycle))

fwrite(metrics, 'project/metrics_18jan2022.csv')