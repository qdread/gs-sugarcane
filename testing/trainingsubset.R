PD <- matrix(NA,nrow=432,ncol=2)
k <- 5
training_size = 0.6

N <- nrow(PD)
fold_ids <- sort(rep_len(1:k, N))
fold_n <- table(fold_ids)
train_use <- rep(TRUE, N)
if (training_size < 1) {
  n_train_remove <- floor((1 - training_size) * fold_n)
  n_train_keep <- fold_n - n_train_remove
  train_use <- rep(rep(c(TRUE, FALSE), k), c(rbind(n_train_keep, n_train_remove)))
}
random_order <- sample(1:N)
fold_ids <- fold_ids[random_order]
train_use <- train_use[random_order]

matrix_ranks <- map_int(1:k, ~ rankMatrix(GD[!fold_ids %in% . & train_use, ]))
training_set_sizes <- ceiling(training_size * (N-fold_n))
rank_deficient <- matrix_ranks < training_set_sizes