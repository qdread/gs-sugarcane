# Testing the ADE model with different implementations of how the fitted values are gotten.
# See if there was an issue.


source('GS_sugarcane_2022_loaddata.R')

#GD, PD, crop_cycle_to_use, trait_to_use, k, marker_density, training_size
GD=geno_mat
PD=pheno_blups
crop_cycle_to_use='PlantCane'
trait_to_use='diam'
k=5
marker_density=1
training_size=1

PD <- PD[trait == trait_to_use, c('Clone', crop_cycle_to_use), with = FALSE]

rank_deficient <- rep(TRUE, k)
while (any(rank_deficient)) {
  N <- nrow(PD)
  fold_ids <- sort(rep_len(1:k, N))
  fold_n <- table(fold_ids)
  
  n_train_remove <- floor((1 - training_size) * fold_n)
  n_train_keep <- fold_n - n_train_remove
  train_use <- rep(rep(c(TRUE, FALSE), k), c(rbind(n_train_keep, n_train_remove)))
  
  random_order <- sample(1:N)
  fold_ids <- fold_ids[random_order]
  train_use <- train_use[random_order]
  
  matrix_ranks <- map_int(1:k, ~ rankMatrix(GD[!fold_ids %in% . & train_use, ]))
  training_set_sizes <- as.integer(map_dbl(1:k, ~ sum(n_train_keep[-.])))
  rank_deficient <- matrix_ranks < training_set_sizes
}

pred_values <- list()

for (fold in 1:k) {
  message('------\nFold ', fold, '\n------')
  # Create training and test set for GD and PD
  train_clones <- PD$Clone[!fold_ids %in% fold & train_use]
  test_clones <- PD$Clone[fold_ids %in% fold]
  
  PD_train <- PD[Clone %in% train_clones]
  PD_test <- PD[Clone %in% test_clones]
  
  GD_train <- GD[match(PD_train$Clone, dimnames(GD)[[1]]), ]
  GD_test <- GD[match(PD_test$Clone, dimnames(GD)[[1]]), ]
  
  Y_train <- PD_train[[2]]
  Y_test <- PD_test[[2]]
  
  
  Y_pred_ADEold <- gs_ADE(Y_train, Y_test, GD_train, GD_test)
 
  
  # Store results in data frame with fold ID, observed phenotype, and 1 column for each model's prediction
  pred_values[[fold]] <- data.frame(fold = fold, Clone = PD_test[['Clone']], Y_obs = Y_test, ADEold = Y_pred_ADEold, ADEnew=Y_pred_ADEnew)
}


ETA_ADE <- list(add=list(Z=Z,K=A), dom=list(Z=Z,K=D), epi=list(Z=Z,K=E))

###############

PD_comb = rbind(PD_train, PD_test)
PD_comb[idx_test, PlantCane := NA]
PD_comb[, CloneD := Clone]
PD_comb[, CloneE := Clone]

A <- A.mat(GD_comb) # additive relationship matrix 
D <- D.mat(GD_comb) # dominance
E <- E.mat(GD_comb) # epistasis

fit_ADE_better <- mmer(fixed = PlantCane ~ 1, random = ~ vs(Clone, Gu = A) + vs(CloneD, Gu = D) + vs(CloneE, Gu = E), rcov = ~ units, data = PD_comb)

U <- with(fit_ADE_better$U, cbind(`u:Clone`[[1]], `u:CloneD`[[1]], `u:CloneE`[[1]]))
B <- fit_ADE_better$Beta$Estimate
Y_pred_all <- rowSums(U) + B

U <- do.call(cbind, fit_ADE_better$U[[1]])
B <- fit_ADE_better$Beta$Estimate
Y_pred_all <- sweep(U, 2, B, `+`)

fit_ADE_better <- mmer(fixed = PlantCane ~ 1, random ~ vs(Clone, Gu = A), rcov = ~ units, data = PD_comb)



#mmer(fixed = model_formula, random = ~ vs(Clone, Gu = A), rcov = ~ units, data = PD_masked)