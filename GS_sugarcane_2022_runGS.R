# GENOMIC SELECTION SUGARCANE: SCRIPT TO RUN GS MODELS
# Author: Quentin Read
# Date: 12 January 2022
# ----------------------------------------------------

source('GS_sugarcane_2022_loaddata.R')

# Run GS models -----------------------------------------------------------

n_iter <- 25 
n_folds <- 5

# Repeat for each:
# - crop cycle (plant cane, 1st ratoon, 2nd ratoon)
# - trait
# - model

# Then in each case do n_iter iterations, and in each iteration the n_folds CV folds.

# For now, just do iterations, models, crop cycles, and traits. Do not vary training size or marker density.

combos <- CJ(iter = 1:n_iter, trait = c(physical_traits, economic_traits), crop_cycle = c('PlantCane', 'Ratoon1', 'Ratoon2'))

# Do the GS. Write observed and predicted phenotypes and prediction accuracy metrics with each iteration.
# If the prediction metric file already exists, skip that iteration (this allows the script to be rerun)
gs_pred_metrics <- future_pmap(combos, function(iter, trait, crop_cycle) {
  metric_file_name <- glue('project/output/metrics_{trait}_{crop_cycle}_{iter}.csv')
  if (!file.exists(metric_file_name)) {
    pred_vals <- gs_all(GD = geno_mat, PD = pheno_blups, 
                        crop_cycle_to_use = crop_cycle, trait_to_use = trait, k = n_folds, marker_density = 1)
    fwrite(pred_vals, glue('project/output/phenotypes_{trait}_{crop_cycle}_{iter}.csv'))
    pred_metrics <- pred_vals[, calc_metrics(Y_obs, Y_pred), by = model]
    fwrite(pred_metrics, metric_file_name)
  } else {
    pred_metrics <- fread(metric_file_name)
  }
  return(pred_metrics)
}, .options = furrr_options(seed = 919))

combos[, metrics := gs_pred_metrics]
results <- unnest_dt(combos, col = metrics, id = .(iter, trait, crop_cycle))

fwrite(results, 'project/output/all_metrics.csv')
