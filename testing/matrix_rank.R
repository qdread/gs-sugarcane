# Test whether we can get rank deficient matrices

N <- 1000
k <- 5
set.seed(77)
ranx <- list()
fullranx <- list()

for (i in 1:N) {
  if (i%%10 == 0) message(i)
  fold_ids <- sample(rep_len(1:k, nrow(PD)))
  ranx[[i]] <- map_int(1:k, ~ rankMatrix(GD[fold_ids == ., ]))
  fullranx[[i]] <- table(fold_ids)
  if (any(ranx[[i]] < fullranx[[i]])) break()
}