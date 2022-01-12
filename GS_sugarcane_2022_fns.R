# Function to split dataset into training and test sets
# Given n (number of individuals in dataset), k (number of folds), and p (proportion of dataset used for training)
split_train_test <- function(k, p, n) {
  train_size <- ceiling(p*n)
  replicate(k, sample(1:n, train_size, replace = FALSE))
}


# Master GS function ------------------------------------------------------

# GS function. This will do one repetition for one trait and a given value of marker density
# Also supply an output directory to write the intermediate files
gs_all <- function(GD, PD, crop_cycle_to_use, trait, k, marker_density, rand_seed, output_dir = NULL) {
  # Get clone ID and trait value for the given crop cycle.
  PD <- PD[crop_cycle == crop_cycle_to_use, c('Clone', trait), with = FALSE]
  
  # If no random seed is supplied, use the current time.
  if (missing(rand_seed)) rand_seed <- floor(as.numeric(Sys.time())%%123456)
  set.seed(rand_seed)
  
  # Subsample markers to given density
  if (marker_density < 1) {
    n_markers <- floor(ncol(GD) * marker_density)
    GD <- GD[, sample(ncol(GD), n_markers, replace = FALSE)]
  }
  
  # Assign individuals to cross-validation folds
  fold_ids <- sample(rep_len(1:k, nrow(PD)))
  
  res <- list()
  
  for (fold in 1:k) {
    message('------\nFold ', fold, '\n------')
    # Create training and test set for GD and PD
    train_clones <- PD$Clone[!fold_ids %in% fold]
    test_clones <- PD$Clone[fold_ids %in% fold]
    
    PD_train <- PD[Clone %in% train_clones]
    PD_test <- PD[!Clone %in% train_clones]
    
    GD_train <- GD[match(PD_train$Clone, dimnames(GD)[[1]]), ]
    GD_test <- GD[match(PD_test$Clone, dimnames(GD)[[1]]), ]
    
    Y_train <- PD_train[[trait]]
    Y_test <- PD_test[[trait]]
    
    # Apply the five models and return the diagnostic values
    message('Running rrBLUP (1 of 5)')
    diags_rrBLUP <- gs_rrBLUP(Y_train, Y_test, GD_train, GD_test)
    message('Running ADE (2 of 5)')
    diags_ADE <- gs_ADE(Y_train, Y_test, GD_train, GD_test)
    message('Running BGLR (3 of 5)')
    diags_BGLR <- gs_BGLR(Y_train, Y_test, GD_train, GD_test)
    message('Running SVM (4 of 5)')
    diags_SVM <- gs_SVM(Y_train, Y_test, GD_train, GD_test)
    message('Running RF (5 of 5)')
    diags_RF <- gs_RF(Y_train, Y_test, GD_train, GD_test)
    
    # Store results in data frame with fold ID and model type
    res[[fold]] <- data.frame(model = c('RRBLUP', 'ADE', 'BGLR', 'SVM', 'RF'),
                              fold = fold,
                              do.call(rbind, list(diags_rrBLUP, diags_ADE, diags_BGLR, diags_SVM, diags_RF)))
  }
 
  # Combine list of data frames together and return
  return(do.call(rbind, res))

}


# Function to calculate prediction accuracy metrics -----------------------

# This function will take the observed and predicted trait (Y) values from any model and return the full set of prediction accuracy metrics
calc_metrics <- function(Y_obs, Y_pred) {
  r <- cor(Y_obs, Y_pred, use = 'complete')
  mod <- lm(Y_obs ~ Y_pred)
  ci <- CI(Y_obs, Y_pred, s = 0.2, top = TRUE)
  rmse <- sqrt(mean((Y_obs - Y_pred)^2))
  
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
  predYr.Ro = predYv.Ro[,1] + rrMod.Ro$beta
  
  # FIXME Here, write the observed and predicted phenotypic values to a file if needed
  calc_metrics(Y_obs = Y_test, Y_pred = predYr.Ro)
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
  A1.2 <- A.mat(GD_comb) # additive relationship matrix 
  D1.2 <- D.mat(GD_comb) # dominance
  E1.2 <- E.mat(GD_comb) # epistasis
  
  Za.2 <- diag(length(Y_comb)) 
  Zd.2 <- diag(length(Y_comb))
  Ze.2 <- diag(length(Y_comb))
  
  rownames(A1.2) <- 1:nrow(A1.2)
  rownames(D1.2) <- 1:nrow(D1.2)
  rownames(E1.2) <- 1:nrow(E1.2)
  
  ETA.AD.2 <- list(add=list(Z=Za.2,K=A1.2), dom=list(Z=Zd.2,K=D1.2), epi=list(Z=Ze.2,K=E1.2))
  ans.AD.2 <- MEMMA(Y=Y_comb, ZETA=ETA.AD.2) 
  
  # FIXME Here write observed and predicted phenotypes to a file if needed
  calc_metrics(Y_obs = Y_test, Y_pred = ans.AD.2$fitted.y[idx_test])
}


# BGLR GS function --------------------------------------------------------

gs_BGLR <- function(Y_train, Y_test, GD_train, GD_test) {

  # Combine train and test sets into single vector (Y) and matrix (GD)
  Y_comb <- c(Y_train, rep(NA, length(Y_test)))
  idx_train <- 1:length(Y_train)
  idx_test <- (1:length(Y_test)) + length(Y_train)
  
  GD_comb <- rbind(GD_train, GD_test)

  M.2 <-tcrossprod(GD_comb)/ncol(GD_comb)

  ETA.2 <-list(list(K = M.2,model = 'RKHS')) 
  fm.RK.2 <- BGLR(y=Y_comb, ETA=ETA.2, response_type="gaussian", nIter=12000, burnIn=2000, verbose = FALSE)
  
  #FIXME Write observed and predicted phenotypes to file if needed
  calc_metrics(Y_obs = Y_test, Y_pred = fm.RK.2$yHat[idx_test])
}


# SVM GS function ---------------------------------------------------------

gs_SVM <- function(Y_train, Y_test, GD_train, GD_test) {
  #FIXME Any tuning parameters?
  svm.model <- svm(GD_train, Y_train, type = 'eps-regression', kernel = 'radial', scale = FALSE)
  svm_pred <- predict(svm.model, GD_test, decision.values = TRUE)
  
  #FIXME Write observed and predicted phenotypes to file if needed
  calc_metrics(Y_obs = Y_test, Y_pred = svm_pred)
}


# Randomforest GS function ------------------------------------------------

gs_RF <- function(Y_train, Y_test, GD_train, GD_test) {
  
  # Adjust values to -1, 0, 1
  GD_train <- GD_train - 1
  GD_test <- GD_test - 1
  
  # FIXME Grid search tuning parameters?
  rf <- randomForest(GD_train, Y_train, tree=5000, mtry=floor(sqrt(ncol(GD_train))), 
                     replace=TRUE, classwt=NULL, nodesize = 1, importance=FALSE, localImp=FALSE, nPerm=1,norm.votes=TRUE, do.trace=FALSE)
  
  rf_pred = predict(rf, newdata=GD_test)
  
  #FIXME Write observed and predicted phenotypes to file if needed
  calc_metrics(Y_obs = Y_test, Y_pred = rf_pred)
}
