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

### test svm tuning

svm_untuned <- svm(GD_train, Y_train, type = 'eps-regression', kernel = 'radial', scale = FALSE)
svmlinear_untuned <- svm(GD_train, Y_train, type = 'eps-regression', kernel = 'linear', scale = FALSE)
svmpoly_untuned <- svm(GD_train, Y_train, type = 'eps-regression', kernel = 'polynomial', scale = FALSE)
svmsig_untuned <- svm(GD_train, Y_train, type = 'eps-regression', kernel = 'sigmoid', scale = FALSE)

gamma1 <- 1/ncol(GD_train)
gamma_grid <- 10^(-2:2) * gamma1
cost_grid <- c(.1, 1, 10)

svm_tuned <- tune(svm, train.x = GD_train, train.y = Y_train, validation.x = GD_test, validation.y = Y_test,
     ranges = list(cost = cost_grid, gamma = gamma_grid),
     type = 'eps-regression', kernel = 'radial', scale = FALSE,
     tunecontrol = tune.control(sampling = 'fix'))

svmlinear_tuned <- tune(svm, train.x = GD_train, train.y = Y_train, validation.x = GD_test, validation.y = Y_test,
                  ranges = list(cost = cost_grid, gamma = gamma_grid),
                  type = 'eps-regression', kernel = 'linear', scale = FALSE,
                  tunecontrol = tune.control(sampling = 'fix'))


Y_pred_untuned <- predict(svm_untuned, GD_test, decision.values = TRUE)
Y_pred_untunedlinear <- predict(svmlinear_untuned, GD_test, decision.values = TRUE)
Y_pred_untunedpoly <- predict(svmpoly_untuned, GD_test, decision.values = TRUE)
Y_pred_untunedsig <- predict(svmsig_untuned, GD_test, decision.values = TRUE)
Y_pred_tuned <- predict(svm_tuned$best.model, GD_test, decision.values = TRUE)
Y_pred_tunedlinear <- predict(svmlinear_tuned$best.model, GD_test, decision.values = TRUE)

calc_metrics(Y_test, Y_pred_untuned)
calc_metrics(Y_test, Y_pred_untunedlinear)
calc_metrics(Y_test, Y_pred_untunedpoly)
calc_metrics(Y_test, Y_pred_untunedsig)
calc_metrics(Y_test, Y_pred_tuned)
calc_metrics(Y_test, Y_pred_tunedlinear)

### test rf tuning

# the randomForest function does not accept NA in response so remove those rows
na_idx <- is.na(Y_train)
Y_train_rf <- Y_train[!na_idx]

# Adjust values to -1, 0, 1
GD_train_rf <- GD_train[!na_idx, ] - 1
GD_test_rf <- GD_test - 1

rf_untuned <- randomForest(GD_train_rf, Y_train_rf, ntree=1000, mtry=floor(sqrt(ncol(GD_train))), 
                           replace=TRUE, classwt=NULL, nodesize = 1, importance=FALSE, localImp=FALSE, nPerm=1, norm.votes=TRUE, do.trace=FALSE)



mtry1 <- floor(sqrt(ncol(GD_train)))
mtry_grid <- c(10, 20, 50, 101, 200)
nodesize_grid <- c(1, 5) # Too low of numbers may be very slow.

rf_tuned <- tune(randomForest, train.x = GD_train_rf, train.y = Y_train_rf, validation.x = GD_test_rf, validation.y = Y_test,
                 ranges = list(mtry = mtry_grid, nodesize = nodesize_grid),
                 ntree = 100, replace=TRUE, classwt=NULL, importance=FALSE, localImp=FALSE, nPerm=1, norm.votes=TRUE, do.trace=TRUE,
                 tunecontrol = tune.control(sampling = 'fix'))

Y_pred_untuned = predict(rf_untuned, newdata=GD_test_rf)
Y_pred_tuned = predict(rf_tuned$best.model, newdata=GD_test_rf)

calc_metrics(Y_test, Y_pred_untuned)
calc_metrics(Y_test, Y_pred_tuned)
