# GENOMIC SELECTION SUGARCANE: SCRIPT TO RUN GS MODELS AT DIFFERENT TRAINING SIZES USING RSLURM
# Author: Quentin Read
# Date: 1 March 2022
# ----------------------------------------------------

source('GS_sugarcane_2022_loaddata_rslurm.R')

# Run GS models -----------------------------------------------------------

n_iter <- 25
n_folds <- 5

# Repeat for each:
# - crop cycle (plant cane, 1st ratoon, 2nd ratoon)
# - trait
# - model
# - training size

# Then in each case do n_iter iterations, and in each iteration the n_folds CV folds.

combos <- CJ(iter = 1:n_iter, trait = c(physical_traits, economic_traits), crop_cycle = c('PlantCane', 'Ratoon1', 'Ratoon2'),
             train_size = c(0.20, 0.30, 0.50, 0.60, 0.80, 0.90))
set.seed(45678)
combos$seed <- round(runif(nrow(combos)) * 76545)

# Do the GS. Write observed and predicted phenotypes and prediction accuracy metrics with each iteration.
# If the prediction metric file already exists, skip that iteration (this allows the script to be rerun)
gs_ts_fun <- function(iter, trait, crop_cycle, train_size, seed) {
  metric_file_name <- glue('/project/qdr/gs_sugarcane/output/TS/metrics_{trait}_{crop_cycle}_TS{train_size}_{iter}.csv')
  if (!file.exists(metric_file_name)) {
    source('/home/quentin.read/GitHub/gs-sugarcane/GS_sugarcane_2022_fns.R')
    set.seed(seed)
    pred_vals <- gs_all(GD = geno_mat, PD = pheno_blups, 
                        crop_cycle_to_use = crop_cycle, trait_to_use = trait, k = n_folds, 
                        marker_density = 1, training_size = train_size, temp_dir = '/90daydata/shared/qdr/gs_sugarcane/temp')
    fwrite(pred_vals, glue('/project/qdr/gs_sugarcane/output/TS/phenotypes_{trait}_{crop_cycle}_TS{train_size}_{iter}.csv'))
    pred_metrics <- pred_vals[, calc_metrics(Y_obs, Y_pred), by = model]
    fwrite(pred_metrics, metric_file_name)
  } else {
    pred_metrics <- fread(metric_file_name)
  }
  return(pred_metrics)
}
tsjob <- slurm_apply(f = gs_ts_fun, params = combos, jobname = 'gsts_rslurm', nodes = 6, cpus_per_node = 20,
                    global_objects = c('geno_mat', 'pheno_blups', 'n_folds'),
                    slurm_options = list(partition = 'medium', time = '7-00:00:00', mem = '80gb'))

# combos[, metrics := gs_pred_metrics]
# results <- unnest_dt(combos, col = metrics, id = .(iter, trait, crop_cycle, train_size))
# 
# fwrite(results, 'project/output/TS_all_metrics.csv')
