# Function to split dataset into training and test sets
# Given n (number of individuals in dataset), k (number of folds), and p (proportion of dataset used for training)
split_train_test <- function(k, p, n) {
  train_size <- ceiling(p*n)
  replicate(k, sample(1:n, train_size, replace = FALSE))
}


# Master GS function ------------------------------------------------------

# GS function. This will do one repetition for one trait and a given value of marker density
# Also supply an output directory to write the intermediate files
gs_all <- function(GD, PD, crop_cycle_to_use, trait, k, marker_density, output_dir) {
  # Get clone ID and trait value for the given crop cycle.
  PD <- PD[crop_cycle == crop_cycle_to_use, c('Clone', trait), with = FALSE]
  
  # If no random seed is supplied, use the current time.
  if (missing(rand_seed)) rand_seed <- floor(as.numeric(Sys.time())%%123456)
  set.seed(rand_seed)
  
  # Subsample markers to given density
  n_markers <- floor(nrow(geno_mat) * marker_density)
  GD <- GD[sample(nrow(GD), n_markers, replace = FALSE), ]
  
  # Assign individuals to cross-validation folds
  fold_ids <- sample(rep_len(1:k, nrow(PD)))
  
  res <- list()
  
  for (fold in 1:k) {
    # Create training and test set for GD and PD
    train_clones <- PD$Clone[!fold_ids %in% fold]
    test_clones <- PD$Clone[fold_ids %in% fold]
    
    PD_train <- PD[Clone %in% train_clones]
    PD_test <- PD[!Clone %in% train_clones]
    
    GD_train <- GD[, match(PD_train$Clone, dimnames(GD)[[2]])]
    GD_test <- GD[, match(PD_test$Clone, dimnames(GD)[[2]])]
    
    Y_train <- PD_train[[trait]]
    Y_test <- PD_test[[trait]]
    
    # Apply the five models and return the diagnostic values
    diags_rrBLUP <- gs_rrBLUP(Y_train, Y_test, GD_train, GD_test)
    diags_ADE <- gs_ADE(Y_train, Y_test, GD_train, GD_test)
    diags_BGLR <- gs_BGLR(Y_train, Y_test, GD_train, GD_test)
    diags_SVM <- gs_SVM(Y_train, Y_test, GD_train, GD_test)
    diags_RF <- gs_RF(Y_train, Y_test, GD_train, GD_test)
    
    # Store results in nested list structure to be combined later
    res[[fold]] <- list(rrBLUP = diags_rrBLUP, ADE = diags_ADE, BGLR = diags_BGLR, SVM = diags_SVM, RF = diags_RF)
  }
 
  # Return nested list 
  return(res)

}


# Function to calculate prediction accuracy metrics -----------------------

# This function will take the observed and predicted trait (Y) values from any model and return the full set of prediction accuracy metrics
calc_metrics <- function(Y_obs, Y_pred) {
  r <- cor(Y_obs, Y_pred, use = 'complete')
  mod <- lm(Y_obs ~ Y_pred)
  ci <- CI(Y_obs, Y_pred, s = 0.2, top = TRUE)
  rmse <- sqrt(mean((Y_obs - Y_pred)^2))
  
  return(list(r = r, intercept = mod$coefficients[1], slope = mod$coefficients[2], CI = ci, RMSE = rmse))
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
  
}


# BGLR GS function --------------------------------------------------------

gs_BGLR <- function(Y_train, Y_test, GD_train, GD_test) {
  # FIXME NEED TO EDIT THIS SO IT WILL WORK

  X.2 <- rbind(Geno_trainingset.trn.size, Geno_testset.trn.size)
  M.2 <-tcrossprod(X.2)/ncol(X.2)
  #y.trn.2[valid.comb.2] <- NA
  
  ETA.2 <-list(list(K=M.2,model='RKHS')) 
  fm.RK.2<-BGLR(y=y.trn.2,ETA=ETA.2,response_type="gaussian", nIter=12000, burnIn=2000)
  
  #FIXME Write observed and predicted phenotypes to file if needed
  calc_metrics(Y_obs = Y_test, Y_pred = fm.RK.2$yHat[valid.comb.2])
  r.gy.RK[t,f] <- cor(fm.RK.2$yHat[valid.comb.2], Y_valid, use="complete")
}


# SVM GS function ---------------------------------------------------------

gs_SVM <- function(Y_train, Y_test, GD_train, GD_test) {
  #FIXME Any tuning parameters?
  svm.model <- svm(GD_train, GD_test, type = 'eps-regression', kernel = 'radial', scale = FALSE)
  svm_pred <- predict(svm.model, GD_test, decision.values = TRUE)
  
  #FIXME Write observed and predicted phenotypes to file if needed
  calc_metrics(Y_obs = Y_test, Y_pred = svm_pred)
}


# Randomforest GS function ------------------------------------------------

gs_RF <- function(Y_train, Y_test, GD_train, GD_test) {
  # FIXME Grid search tuning parameters?
  rf <- randomForest(GD_train, Y_train, tree=5000, mtry=floor(sqrt(ncol(GD_train))), 
                     replace=TRUE, classwt=NULL, nodesize = 1, importance=FALSE, localImp=FALSE, nPerm=1,norm.votes=TRUE, do.trace=FALSE)
  
  pred.rf = predict(rf, newdata=Geno_testset.rrblup)
}
