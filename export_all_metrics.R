# GENOMIC SELECTION SUGARCANE: Combine results and write to single file to export (for all four analyses)
# Author: Quentin Read
# Date: 03 March 2022
# -------------------------------------------------------------------------------------------------------

library(data.table)
library(glue)

physical_traits <- c("stkwt_kg", "diam", "Brix", "Fiber", "Pol", "Sucrose", "Purity", "stalk_ha")
economic_traits <- c("TCH", "TRS", "CRS", "TSH", "EI")

### GS

combos_gs <- CJ(iter = 1:25, trait = c(physical_traits, economic_traits), crop_cycle = c('PlantCane', 'Ratoon1', 'Ratoon2'))

combos_gs[, metric_file_name := with(combos_gs, glue('project/output/GS/metrics_{trait}_{crop_cycle}_{iter}.csv'))]
metrics_gs <- combos_gs[, fread(metric_file_name), by = .(iter, trait, crop_cycle)]

fwrite(metrics_gs, 'project/output/GS_all_metrics.csv')

combos_gs[, phenotype_file_name := with(combos_gs, glue('project/output/GS/phenotypes_{trait}_{crop_cycle}_{iter}.csv'))]
phenotypes_gs <- combos_gs[, fread(phenotype_file_name), by = .(iter, trait, crop_cycle)]
phenotypes_gs <- phenotypes_gs[, .(trait, crop_cycle, iter, model, fold, Clone, Y_obs, Y_pred)]

fwrite(phenotypes_gs, 'project/output/GS_observed_predicted_phenotypes.csv')

### MD

combos_md <- CJ(iter = 1:25, trait = c(physical_traits, economic_traits), crop_cycle = c('PlantCane', 'Ratoon1', 'Ratoon2'), 
             marker_density = c(0.20, 0.30, 0.50, 0.60, 0.80, 0.90))

combos_md[, metric_file_name := with(combos_md, glue('project/output/MD/metrics_{trait}_{crop_cycle}_MD{marker_density}_{iter}.csv'))]
metrics_md <- combos_md[, fread(metric_file_name), by = .(iter, trait, crop_cycle, marker_density)]

fwrite(metrics_md, 'project/output/MD_all_metrics.csv')

combos_md[, phenotype_file_name := with(combos_md, glue('project/output/MD/phenotypes_{trait}_{crop_cycle}_MD{marker_density}_{iter}.csv'))]
phenotypes_md <- combos_md[, fread(phenotype_file_name), by = .(iter, trait, crop_cycle, marker_density)]
phenotypes_md <- phenotypes_md[, .(trait, crop_cycle, marker_density, iter, model, fold, Clone, Y_obs, Y_pred)]

fwrite(phenotypes_md, 'project/output/MD_observed_predicted_phenotypes.csv')

### TS

combos_ts <- CJ(iter = 1:25, trait = c(physical_traits, economic_traits), crop_cycle = c('PlantCane', 'Ratoon1', 'Ratoon2'), 
                training_size = c(0.20, 0.30, 0.50, 0.60, 0.80, 0.90))

combos_ts[, metric_file_name := with(combos_ts, glue('project/output/TS/metrics_{trait}_{crop_cycle}_TS{training_size}_{iter}.csv'))]
metrics_ts <- combos_ts[, fread(metric_file_name), by = .(iter, trait, crop_cycle, training_size)]

fwrite(metrics_ts, 'project/output/TS_all_metrics.csv')

combos_ts[, phenotype_file_name := with(combos_ts, glue('project/output/TS/phenotypes_{trait}_{crop_cycle}_TS{training_size}_{iter}.csv'))]
phenotypes_ts <- combos_ts[, fread(phenotype_file_name), by = .(iter, trait, crop_cycle, training_size)]
phenotypes_ts <- phenotypes_ts[, .(trait, crop_cycle, training_size, iter, model, fold, Clone, Y_obs, Y_pred)]

fwrite(phenotypes_ts, 'project/output/TS_observed_predicted_phenotypes.csv')

### TAGS

combos_tags <- CJ(iter = 1:25, trait = c(physical_traits, economic_traits), crop_cycle = c('Ratoon1', 'Ratoon2'))

# Metrics are already written out.

combos_tags[, phenotype_file_name := with(combos_tags, glue('project/output/TAGS/phenotypes_{trait}_{crop_cycle}_{iter}.csv'))]
phenotypes_tags <- combos_tags[, fread(phenotype_file_name), by = .(iter, trait, crop_cycle)]
phenotypes_tags <- phenotypes_tags[, .(trait, crop_cycle, iter, fold, Clone, Y_obs, Y_pred)]

fwrite(phenotypes_tags, 'project/output/TAGS_observed_predicted_phenotypes.csv')
