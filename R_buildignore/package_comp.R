source("spEigen.R")

#--------------------#
# Libraries required #
library(mvtnorm) # rmvnorm function for data generation
library('elasticnet')

set.seed(42)

#------------#
# Parameters #
m <- 300 # dimension
n <- 200 # number of samples
q <- 3 # number of sparse eigenvectors to be estimated
sp_card <- 0.1*m # cardinality of the sparse eigenvectors
rho <- 0.5

#-----------------#
# True Covariance #
# Sparse eigenvectors
V <- matrix(rnorm(m^2), ncol = m)
tmp <- matrix(0, m, q)

for (i in 1:max(q, 2)) {
  ind1 <- (i-1)*sp_card + 1
  ind2 <- i*sp_card
  tmp[ind1:ind2, i] = 1/sqrt(sp_card)
  V[, i] <- tmp[, i]
}

V <- qr.Q(qr(V)) # orthogonalization, but keep the first eigenvectors to be same as V

# Eigenvalues
vl <- rep(1, m)
for (i in 1:q) {
  vl[i] <- 100*(q + 1 - i)
}

# Covariance matrix
R <- V %*% diag(vl) %*% t(V)

#-------------#
# Data Matrix #
X <- rmvnorm(n = n, mean = rep(0, m), sigma = R) # random data with underlying sparse structure
X <- X - matrix(rep(colMeans(X), n), nrow = n, byrow = T) # center the data

#-------------------------------#
# Sparse Eigenvector Extraction #
ptm <- proc.time()
res_spca <- spca(cov(X),K=q,type="Gram",sparse="penalty",trace=TRUE,para=c(.4,.4,.4))
time_spca <- (proc.time() - ptm)[3]

ptm <- proc.time()
res_sparse <- spEigen(X, q, rho, data = TRUE)
time_sparse <- (proc.time() - ptm)[3]

#------#
# Time #
cat('SPCA running time: ', time_spca, '\n')
cat('spEigen running time: ', time_sparse, '\n')

#----------#
# Recovery #
recovery_spca <- abs(diag(Conj(t(res_spca$loadings)) %*% V[, 1:q])) # recovery
print(recovery_spca)
recovery_sparse <- abs(diag(Conj(t(res_sparse$vectors)) %*% V[, 1:q])) # recovery
print(recovery_sparse)


#-------#
# Plots #
par(mfcol = c(3, 2))
plot(res_sparse$vectors[, 1]*sign(res_sparse$vectors[1, 1]), main = "First Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(V[, 1]*sign(V[1, 1]), col = "red")
plot(res_sparse$vectors[, 2]*sign(res_sparse$vectors[sp_card+1, 2]), main = "Second Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
plot(res_sparse$vectors[, 3]*sign(res_sparse$vectors[2*sp_card+1, 3]), main = "Third Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")

plot(res_spca$loadings[, 1]*sign(res_spca$loadings[1, 1]), main = "First Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(V[, 1]*sign(V[1, 1]), col = "red")
plot(res_spca$loadings[, 2]*sign(res_spca$loadings[sp_card+1, 2]), main = "Second Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
plot(res_spca$loadings[, 3]*sign(res_spca$loadings[2*sp_card+1, 3]), main = "Third Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")
