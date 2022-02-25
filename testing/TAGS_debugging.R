# TAGS debugging....

combos_short <- combos[iter==1]

for (i in 1:nrow(combos_short)) {
  
  message("Row ", i)  
  
  pred_vals <- trait_assisted_gs(
    GD = geno_mat, PD = pheno_blups, crop_cycle_to_use = combos_short$crop_cycle[i], trait_to_use = combos_short$trait[i], k = n_folds
  )
  
}