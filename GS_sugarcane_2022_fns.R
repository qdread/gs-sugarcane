# GENOMIC SELECTION SUGARCANE: SOURCE FILE FOR FUNCTIONS
# Author: Quentin Read
# Date: 12 January 2022
# ------------------------------------------------------


# Function to find BLUPs by trait -----------------------------------------

# Crop cycle is fit as a fixed effect.
# If a trait was only measured during one crop cycle, do not fit crop cycle in the model.
blup_trait <- function(dat) {
  if (length(unique(dat$crop_cycle)) > 1) {
    lmm <- lmer(value ~ 0 + crop_cycle + (1|Clone) + (1|Row) + (1|Column), data = dat,
                control = lmerControl(optimizer = 'bobyqa'))
    blup <- outer(ranef(lmm)[['Clone']][['(Intercept)']], fixef(lmm), `+`)
    data.frame(Clone = row.names(ranef(lmm)[['Clone']]), blup)
  } else {
    lmm <- lmer(value ~ 1 + (1|Clone) + (1|Row) + (1|Column), data = dat,
                control = lmerControl(optimizer = 'bobyqa'))
    blup <- outer(ranef(lmm)[['Clone']][['(Intercept)']], fixef(lmm), `+`)
    data.frame(Clone = row.names(ranef(lmm)[['Clone']]), blup, as.numeric(NA), as.numeric(NA))
  }
}

# Master GS function ------------------------------------------------------

# This will do one repetition for one trait and a given value of marker density. Return observed and predicted values.
gs_all <- function(GD, PD, crop_cycle_to_use, trait_to_use, k, marker_density) {
  # Get clone ID and trait value for the given crop cycle.
  PD <- PD[trait == trait_to_use, c('Clone', crop_cycle_to_use), with = FALSE]
  
  # Subsample markers to given density
  if (marker_density < 1) {
    n_markers <- floor(ncol(GD) * marker_density)
    GD <- GD[, sample(ncol(GD), n_markers, replace = FALSE)]
  }
  
  # Assign individuals to cross-validation folds
  # Check whether all training sets are full rank matrices to avoid error
  matrix_ranks <- rep(0, k)
  training_set_sizes <- table(fold_ids)
  while (any(matrix_ranks < training_set_sizes)) {
    fold_ids <- sample(rep_len(1:k, nrow(PD)))
    matrix_ranks <- map_int(1:k, ~ rankMatrix(GD[fold_ids == ., ]))
  }
  
  pred_values <- list()
  
  for (fold in 1:k) {
    message('------\nFold ', fold, '\n------')
    # Create training and test set for GD and PD
    train_clones <- PD$Clone[!fold_ids %in% fold]
    test_clones <- PD$Clone[fold_ids %in% fold]
    
    PD_train <- PD[Clone %in% train_clones]
    PD_test <- PD[!Clone %in% train_clones]
    
    GD_train <- GD[match(PD_train$Clone, dimnames(GD)[[1]]), ]
    GD_test <- GD[match(PD_test$Clone, dimnames(GD)[[1]]), ]
    
    Y_train <- PD_train[[2]]
    Y_test <- PD_test[[2]]
    
    # Apply the five models and return the predicted values
    message('Running rrBLUP (1 of 7)')
    Y_pred_rrBLUP <- gs_rrBLUP(Y_train, Y_test, GD_train, GD_test)
    message('Running ADE (2 of 7)')
    Y_pred_ADE <- gs_ADE(Y_train, Y_test, GD_train, GD_test)
    message('Running RKHS (3 of 7)')
    Y_pred_RKHS <- gs_RKHS(Y_train, Y_test, GD_train, GD_test)
    message('Running BayesA (4 of 7)')
    Y_pred_BayesA <- gs_Bayes(Y_train, Y_test, GD_train, GD_test, bayes_model = 'BayesA')
    message('Running BayesB (5 of 7)')
    Y_pred_BayesB <- gs_Bayes(Y_train, Y_test, GD_train, GD_test, bayes_model = 'BayesB')
    message('Running SVM (6 of 7)')
    Y_pred_SVM <- gs_SVM(Y_train, Y_test, GD_train, GD_test)
    message('Running RF (7 of 7)')
    Y_pred_RF <- gs_RF(Y_train, Y_test, GD_train, GD_test)
    
    # Store results in data frame with fold ID, observed phenotype, and 1 column for each model's prediction
    pred_values[[fold]] <- data.frame(fold = fold, Clone = PD_test[['Clone']], Y_obs = Y_test, rrBLUP = Y_pred_rrBLUP, ADE = Y_pred_ADE, RKHS = Y_pred_RKHS, SVM = Y_pred_SVM, RF = Y_pred_RF)
  }
  
  # Combine observed and predicted values for the folds
  pred_values_allfolds <- do.call(rbind, pred_values) |> 
    setDT() |>
    melt(id.vars = c('fold', 'Clone', 'Y_obs'), variable.name = 'model', value.name = 'Y_pred')

  return(pred_values_allfolds)
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

  # Combine train and test sets into single vector (Y) and matrix (GD)
  Y_comb <- c(Y_train, rep(NA, length(Y_test)))
  idx_train <- 1:length(Y_train)
  idx_test <- (1:length(Y_test)) + length(Y_train)
  
  GD_comb <- rbind(GD_train, GD_test) - 1 # adjust to -1, 0, 1 values
  
  # FIXME the argument shrink = TRUE was removed because it is not in the current up to date sommer package. Check if this is OK.
  A <- A.mat(GD_comb) # additive relationship matrix 
  D <- D.mat(GD_comb) # dominance
  E <- E.mat(GD_comb) # epistasis
  
  Z <- diag(length(Y_comb)) 

  rownames(A) <- 1:nrow(A)
  rownames(D) <- 1:nrow(D)
  rownames(E) <- 1:nrow(E)
  
  ETA_ADE <- list(add=list(Z=Z,K=A), dom=list(Z=Z,K=D), epi=list(Z=Z,K=E))
  fit_ADE <- MEMMA(Y=Y_comb, ZETA=ETA_ADE) 
  
  Y_pred = fit_ADE$fitted.y[idx_test]
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
  
  Y_pred = fit_Bayes$yHat[idx_test]
  return(Y_pred)
}

# SVM GS function ---------------------------------------------------------

gs_SVM <- function(Y_train, Y_test, GD_train, GD_test) {
  #FIXME Any tuning parameters?
  svm.model <- svm(GD_train, Y_train, type = 'eps-regression', kernel = 'radial', scale = FALSE)
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
  
  # FIXME Grid search tuning parameters?
  rf <- randomForest(GD_train, Y_train, ntree=5000, mtry=floor(sqrt(ncol(GD_train))), 
                     replace=TRUE, classwt=NULL, nodesize = 1, importance=FALSE, localImp=FALSE, nPerm=1, norm.votes=TRUE, do.trace=FALSE)
  
  Y_pred = predict(rf, newdata=GD_test)
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

