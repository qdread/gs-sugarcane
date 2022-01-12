# Test sommer
data(DT_cpdata)

head(DT_cpdata)
GT_cpdata[1:5,1:5]

Y <- DT_cpdata[, c('color','Yield')]
Yscale <- scale(Y)
Za <- diag(dim(Y)[1])
A <- A.mat(GT_cpdata)
D <- D.mat(GT_cpdata)
E <- E.mat(GT_cpdata)

eta_a <- list(add=list(Z=Za,K=A))
ans_a <- MEMMA(Y=Yscale, ZETA=eta_a)
