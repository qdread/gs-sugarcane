# Test whether tuning parameters matter for the models run in BGLR

### setup
PD = pheno_blups
GD = geno_mat
k = 5
fold = 1

PD <- PD[trait == 'stkwt_kg', c('Clone', 'Ratoon1'), with = FALSE]


# Assign individuals to cross-validation folds
# Check whether all training sets are full rank matrices to avoid error
rank_deficient <- rep(TRUE, k)
while (any(rank_deficient)) {
  fold_ids <- sample(rep_len(1:k, nrow(PD)))
  matrix_ranks <- map_int(1:k, ~ rankMatrix(GD[!fold_ids %in% ., ]))
  training_set_sizes <- as.integer(nrow(PD) - table(fold_ids))
  rank_deficient <- matrix_ranks < training_set_sizes
}

train_clones <- PD$Clone[!fold_ids %in% fold]
test_clones <- PD$Clone[fold_ids %in% fold]

PD_train <- PD[Clone %in% train_clones]
PD_test <- PD[!Clone %in% train_clones]

GD_train <- GD[match(PD_train$Clone, dimnames(GD)[[1]]), ]
GD_test <- GD[match(PD_test$Clone, dimnames(GD)[[1]]), ]

Y_train <- PD_train[[2]]
Y_test <- PD_test[[2]]

### test tuning

# default df0=5, R2=.5, S0 is calculated so that prior mode of resud variance = var(y)*R2

var(Y_train) * 0.5

Y_pred_BayesC <- gs_Bayes(Y_train,Y_test,GD_train,GD_test,bayes_model="BayesC")
calc_metrics(Y_test, Y_pred_BayesC)

Y_pred_BL <- gs_Bayes(Y_train,Y_test,GD_train,GD_test,bayes_model="BL")
calc_metrics(Y_test, Y_pred_BL)

Y_pred_BRR <- gs_Bayes(Y_train,Y_test,GD_train,GD_test,bayes_model="BRR")
calc_metrics(Y_test, Y_pred_BRR)

Y_pred_BayesA <- gs_Bayes(Y_train,Y_test,GD_train,GD_test,bayes_model="BayesA")
calc_metrics(Y_test, Y_pred_BayesA)


# According to the vignette it looks like the df in BayesA and BayesB is really the only thing worth modifying

# Generate unique ID for temp files
UUID <- UUIDgenerate()

# Combine train and test sets into single vector (Y) and matrix (GD)
Y_comb <- c(Y_train, rep(NA, length(Y_test)))
idx_test <- (1:length(Y_test)) + length(Y_train)

GD_comb <- rbind(GD_train, GD_test)

ETA_BayesA <-list(list(X = GD_comb, model = 'BayesA')) 
ETA_BayesB <-list(list(X = GD_comb, model = 'BayesB')) 

metr_A <- list()
metr_B <- list()

for (i in 1:10) {
  message("DF", i)
  fit_BayesA <- BGLR(y = Y_comb, ETA = ETA_BayesA, df0 = i, response_type = "gaussian", nIter = 12000, burnIn = 2000, verbose = TRUE, saveAt = glue('temp/{UUID}'))
  fit_BayesB <- BGLR(y = Y_comb, ETA = ETA_BayesB, df0 = i, response_type = "gaussian", nIter = 12000, burnIn = 2000, verbose = TRUE, saveAt = glue('temp/{UUID}'))
  
  Y_pred_A = fit_BayesA$yHat[idx_test]
  Y_pred_B = fit_BayesB$yHat[idx_test]
  
  metr_A[[i]] <- calc_metrics(Y_test, Y_pred_A)
  metr_B[[i]] <- calc_metrics(Y_test, Y_pred_B)
}

metr_A <- data.frame(df=1:10, do.call(rbind, metr_A))
metr_B <- data.frame(df=1:10, do.call(rbind, metr_B))

with(metr_A, plot(df, RMSE))
with(metr_B, plot(df, RMSE)) # Looks like it doesn't make any difference
