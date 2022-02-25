# GENOMIC SELECTION SUGARCANE: SOURCE FILE FOR FUNCTIONS
# Author: Quentin Read
# Date: 12 January 2022
# ------------------------------------------------------


# Function to find BLUPs by trait -----------------------------------------

# Crop cycle is fit as a fixed effect if crop_cycle_fixef = TRUE.
blup_trait <- function(dat, crop_cycle_fixef = TRUE) {
    if (crop_cycle_fixef) {
      lmm <- lmer(value ~ 0 + crop_cycle + (1|Clone) + (1|Row) + (1|Column), data = dat,
                control = lmerControl(optimizer = 'bobyqa'))
    } else {
      lmm <- lmer(value ~ (1|Clone) + (1|Row) + (1|Column), data = dat,
                  control = lmerControl(optimizer = 'bobyqa'))
    }
    blup <- outer(ranef(lmm)[['Clone']][['(Intercept)']], fixef(lmm), `+`)
    data.frame(Clone = row.names(ranef(lmm)[['Clone']]), blup)
}

# Master GS function ------------------------------------------------------

# This will do one repetition for one trait and a given value of marker density. Return observed and predicted values.
gs_all <- function(GD, PD, crop_cycle_to_use, trait_to_use, k, marker_density, training_size) {
  # Get clone ID and trait value for the given crop cycle.
  PD <- PD[trait == trait_to_use, c('Clone', crop_cycle_to_use), with = FALSE]
  
  # Subsample markers to given density
  if (marker_density < 1) {
    n_markers <- floor(ncol(GD) * marker_density)
    GD <- GD[, sample(ncol(GD), n_markers, replace = FALSE)]
  }
  
  # Assign individuals to cross-validation folds
  # Check whether all training sets are full rank matrices to avoid error
  # Here also subset it further if training size < 1 (assign NA to the ones not in the subset)
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
    
    # Apply the five models and return the predicted values
    message('Running rrBLUP (1 of 9)')
    Y_pred_rrBLUP <- gs_rrBLUP(Y_train, Y_test, GD_train, GD_test)
    message('Running ADE (2 of 9)')
    Y_pred_ADE <- gs_ADE(Y_train, Y_test, GD_train, GD_test)
    message('Running RKHS (3 of 9)')
    Y_pred_RKHS <- gs_RKHS(Y_train, Y_test, GD_train, GD_test)
    message('Running BayesA (4 of 9)')
    Y_pred_BayesA <- gs_Bayes(Y_train, Y_test, GD_train, GD_test, bayes_model = 'BayesA')
    message('Running BayesB (5 of 9)')
    Y_pred_BayesB <- gs_Bayes(Y_train, Y_test, GD_train, GD_test, bayes_model = 'BayesB')
    message('Running SVM radial (6 of 9)')
    Y_pred_SVMradial <- gs_SVM(Y_train, Y_test, GD_train, GD_test, kernel_type = 'radial')
    message('Running SVM linear (7 of 9)')
    Y_pred_SVMlinear <- gs_SVM(Y_train, Y_test, GD_train, GD_test, kernel_type = 'linear')
    message('Running SVM sigmoid (8 of 9)')
    Y_pred_SVMsigmoid <- gs_SVM(Y_train, Y_test, GD_train, GD_test, kernel_type = 'sigmoid')
    message('Running RF (9 of 9)')
    Y_pred_RF <- gs_RF(Y_train, Y_test, GD_train, GD_test)
    
    # Store results in data frame with fold ID, observed phenotype, and 1 column for each model's prediction
    pred_values[[fold]] <- data.frame(fold = fold, Clone = PD_test[['Clone']], Y_obs = Y_test, rrBLUP = Y_pred_rrBLUP, ADE = Y_pred_ADE, RKHS = Y_pred_RKHS, BayesA = Y_pred_BayesA, BayesB = Y_pred_BayesB, SVMradial = Y_pred_SVMradial, SVMlinear = Y_pred_SVMlinear, SVMsigmoid = Y_pred_SVMsigmoid, RF = Y_pred_RF)
  }
  
  # Combine observed and predicted values for the folds
  pred_values_allfolds <- do.call(rbind, pred_values) |> 
    setDT() |>
    melt(id.vars = c('fold', 'Clone', 'Y_obs'), variable.name = 'model', value.name = 'Y_pred')

  return(pred_values_allfolds)
}


# Trait-assisted GS -------------------------------------------------------

trait_assisted_gs <- function(GD, PD, crop_cycle_to_use, trait_to_use, k) {
  # If crop cycle to predict is Ratoon1, use PlantCane to predict
  # If crop cycle to predict is Ratoon2, use PlantCane and Ratoon1 to predict
  PD <- PD[trait == trait_to_use]
  PD <- dcast(PD, Clone ~ crop_cycle, value.var = 'BLUP')
  if (crop_cycle_to_use == 'Ratoon1') model_formula <- as.formula(cbind(PlantCane, Ratoon1) ~ 1)
  if (crop_cycle_to_use == 'Ratoon2') model_formula <- as.formula(cbind(PlantCane, Ratoon1, Ratoon2) ~ 1)
  
  # Assign individuals to cross-validation folds
  fold_ids <- sample(rep_len(1:k, nrow(PD)))

  pred_values <- list()
  
  for (fold in 1:k) {
    message('------\nFold ', fold, '\n------')
    # Create training and test set for PD
    train_clones <- PD$Clone[!fold_ids %in% fold]
    test_clones <- PD$Clone[fold_ids %in% fold]
    
    GD_train <- GD[match(train_clones, dimnames(GD)[[1]]), ]
    GD_test <- GD[match(test_clones, dimnames(GD)[[1]]), ]
    GD_combined <- rbind(GD_train, GD_test)
    
    A <- A.mat(GD_combined)
    
    # Set values for the crop cycle to predict to NA in test set only
    Y_obs <- PD[Clone %in% test_clones][[crop_cycle_to_use]]
    PD_masked <- copy(PD)
    PD_masked[Clone %in% test_clones, (crop_cycle_to_use) := NA]
    
    # Fit model and extract BLUPs (predicted values) for only the crop cycle to predict and only the test set.
    fit_trait_gs <- mmer(fixed = model_formula, random = ~ vs(Clone, Gu = A), rcov = ~ units, data = PD_masked)

    U <- do.call(cbind, fit_trait_gs$U[[1]])
    B <- fit_trait_gs$Beta$Estimate
    Y_pred_all <- sweep(U, 2, B, `+`)
    Y_pred <- Y_pred_all[test_clones, crop_cycle_to_use]
    
    # Store results in data frame with fold ID, observed phenotype, and 1 column for each model's prediction
    pred_values[[fold]] <- data.table(fold = fold, Clone = test_clones, Y_obs = Y_obs, Y_pred = Y_pred)
  }
  
  do.call(rbind, pred_values)
  
}


# Function to calculate prediction accuracy metrics -----------------------

# This function will take the observed and predicted trait (Y) values from any model and return the full set of prediction accuracy metrics
calc_metrics <- function(Y_obs, Y_pred) {
  r <- cor(Y_obs, Y_pred, use = 'complete')
  mod <- lm(Y_obs ~ Y_pred)
  ci <- CI(Y_obs, Y_pred, s = 0.2, top = TRUE)
  rmse <- sqrt(mean((Y_obs - Y_pred)^2, na.rm = TRUE))
  
  return(data.frame(r = r, intercept = mod$coefficients[1], slope = mod$coefficients[2], CI = ci, RMSE = rmse))
}



# CI function -------------------------------------------------------------

# Coincidence index function taken from Marcus' code, and modified to not need the taxon IDs
CI <- function(x, y, s, top = TRUE) {
  
  size <- ceiling(length(x) * s)
  x_order <- order(x, decreasing = top)
  y_order <- order(y, decreasing = top)
  
  both <- sum(x_order[1:size] %in% y_order[1:size])
  random <- both * s
  ci <- (both - random)/(size - random)
  
  return(ci)
}

# rrBLUP GS function ------------------------------------------------------

# Code here is modified from Marcus' code in GS.Marker.Density.... script as well as the Sugarcane.GS.RRBLUP... script

gs_rrBLUP <- function(Y_train, Y_test, GD_train, GD_test) {
  
  GD_train_rrBLUP <- GD_train - 1
  GD_test_rrBLUP <- GD_test - 1

  rrMod.Ro <- mixed.solve(Y_train, X = NULL, Z = GD_train_rrBLUP, K = NULL, SE = FALSE, return.Hinv = FALSE)
  mEff.Ro <- rrMod.Ro$u
  e.Ro = as.matrix(mEff.Ro)
  predYv.Ro = GD_test_rrBLUP %*% e.Ro
  Y_pred = predYv.Ro[,1] + c(rrMod.Ro$beta)
  
  return(Y_pred)
}


# ADE GS function ---------------------------------------------------------

# Again code is modified from Marcus' code.

gs_ADE <- function(Y_train, Y_test, GD_train, GD_test) {

  # Combine train and test sets into single data frame with clone IDs and matrix (GD)
  GD_comb <- rbind(GD_train, GD_test) - 1 # adjust to -1, 0, 1 values
  
  Y_comb <- c(Y_train, rep(NA, length(Y_test)))
  ids <- dimnames(GD_comb)[[1]]
  PD_comb <- data.frame(idA = ids, idD = ids, idE = ids, Y = Y_comb)
  
  idx_train <- 1:length(Y_train)
  idx_test <- (1:length(Y_test)) + length(Y_train)
  
  
  
  # Note: The argument shrink = TRUE was removed because it is not in the current up to date sommer package.
  A <- A.mat(GD_comb) # additive relationship matrix 
  D <- D.mat(GD_comb) # dominance
  E <- E.mat(GD_comb) # epistasis
  
  # Fit model and extract BLUPs (predicted values) for the test set.
  fit_ADE <- mmer(fixed = Y ~ 1, random = ~ vs(idA, Gu = A) + vs(idD, Gu = D) + vs(idE, Gu = E), 
                  rcov = ~ units, data = PD_comb)
  
  U <- with(fit_ADE$U, cbind(`u:idA`[[1]], `u:idD`[[1]], `u:idE`[[1]]))
  B <- fit_ADE$Beta$Estimate
  Y_pred_all <- rowSums(U) + B

  Y_pred <- Y_pred_all[idx_test]
  return(Y_pred)
}


# RKHS GS function --------------------------------------------------------

# RKHS model implemented in BGLR package

gs_RKHS <- function(Y_train, Y_test, GD_train, GD_test) {
  # Generate unique ID for temp files
  UUID <- UUIDgenerate()

  # Combine train and test sets into single vector (Y) and matrix (GD)
  Y_comb <- c(Y_train, rep(NA, length(Y_test)))
  idx_test <- (1:length(Y_test)) + length(Y_train)
  
  GD_comb <- rbind(GD_train, GD_test)

  M <- tcrossprod(GD_comb)/ncol(GD_comb)

  ETA_RKHS <-list(list(K = M, model = 'RKHS')) 
  fit_RKHS <- BGLR(y = Y_comb, ETA = ETA_RKHS, response_type = "gaussian", nIter = 12000, burnIn = 2000, verbose = FALSE, saveAt = glue('temp/{UUID}'))
  
  Y_pred = fit_RKHS$yHat[idx_test]
  return(Y_pred)
}


# Bayes A/B GS function ---------------------------------------------------

# As implemented in BGLR package
# Specify bayes_model 'BayesA' or 'BayesB'
gs_Bayes <- function(Y_train, Y_test, GD_train, GD_test, bayes_model) {
  # Generate unique ID for temp files
  UUID <- UUIDgenerate()
  
  # Combine train and test sets into single vector (Y) and matrix (GD)
  Y_comb <- c(Y_train, rep(NA, length(Y_test)))
  idx_test <- (1:length(Y_test)) + length(Y_train)
  
  GD_comb <- rbind(GD_train, GD_test)
  
  ETA_Bayes <-list(list(X = GD_comb, model = bayes_model)) 
  fit_Bayes <- BGLR(y = Y_comb, ETA = ETA_Bayes, response_type = "gaussian", nIter = 12000, burnIn = 2000, verbose = FALSE, saveAt = glue('temp/{UUID}'))
  
  Y_pred <- fit_Bayes$yHat[idx_test]
  return(Y_pred)
}

# SVM GS function ---------------------------------------------------------

gs_SVM <- function(Y_train, Y_test, GD_train, GD_test, kernel_type) {
  
  svm.model <- svm(GD_train, Y_train, type = 'eps-regression', kernel = kernel_type, scale = FALSE)
  Y_pred <- predict(svm.model, GD_test, decision.values = TRUE)
  
  return(Y_pred)
}


# Randomforest GS function ------------------------------------------------

gs_RF <- function(Y_train, Y_test, GD_train, GD_test) {
  
  # the randomForest function does not accept NA in response so remove those rows
  na_idx <- is.na(Y_train)
  Y_train <- Y_train[!na_idx]
  
  # Adjust values to -1, 0, 1
  GD_train <- GD_train[!na_idx, ] - 1
  GD_test <- GD_test - 1
  
  rf <- randomForest(GD_train, Y_train, ntree = 5000, mtry = floor(sqrt(ncol(GD_train))), nodesize = 1, nPerm = 1,
                     replace = TRUE, classwt = NULL, importance = FALSE, localImp = FALSE, norm.votes = TRUE, do.trace = FALSE)
  
  Y_pred <- predict(rf, newdata = GD_test)
  return(Y_pred)
}


# Utility to unnest data table --------------------------------------------

unnest_dt <- function(dt, col, id, id_vars = NULL) {
  stopifnot(is.data.table(dt))
  
  if (missing(id_vars)) {
    by <-substitute(id)
    col <-substitute(unlist(col, recursive = FALSE))
    
    dt[, eval(col), by = eval(by)]
  } else {
    col <-substitute(unlist(col, recursive = FALSE))
    dt[, eval(col), by = id_vars]
  }
}

# Function to split data into training and test ---------------------------

# (currently not being used)

# Given n (number of individuals in dataset), k (number of folds), and p (proportion of dataset used for training)
split_train_test <- function(k, p, n) {
  train_size <- ceiling(p*n)
  replicate(k, sample(1:n, train_size, replace = FALSE))
}

