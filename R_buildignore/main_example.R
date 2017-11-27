source("spEigen.R")
source("spEigenCov.R")


#--------------------#
# Libraries required #
library(mvtnorm) # rmvnorm function for data generation

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
res_standard <- eigen(cov(X))
res_sparse <- spEigen(X, q, rho, data = TRUE)
if (n > m) {
  res_sparseCov <- spEigenCov(cov(X), q, rho)
}

#-------#
# Plots #
recovery_standard <- abs(diag(Conj(t(res_standard$vectors[, 1:q])) %*% V[, 1:q])) # recovery
print(recovery_standard)
recovery_sparse <- abs(diag(Conj(t(res_sparse$vectors)) %*% V[, 1:q])) # recovery
print(recovery_sparse)
if (n > m) {
  recovery_sparseCov <- abs(diag(Conj(t(res_sparseCov$vectors)) %*% V[, 1:q])) # recovery
  print(recovery_sparseCov)
}

if (n <= m) {
  par(mfcol = c(3, 2))
  plot(res_sparse$vectors[, 1]*sign(res_sparse$vectors[1, 1]), main = "First Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
  lines(V[, 1]*sign(V[1, 1]), col = "red")
  plot(res_sparse$vectors[, 2]*sign(res_sparse$vectors[sp_card+1, 2]), main = "Second Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
  lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
  plot(res_sparse$vectors[, 3]*sign(res_sparse$vectors[2*sp_card+1, 3]), main = "Third Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
  lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")

  plot(res_standard$vectors[, 1]*sign(res_standard$vectors[1, 1]), main = "First Eigenvector", xlab = "Index", ylab = "", type = "h")
  lines(V[, 1]*sign(V[1, 1]), col = "red")
  plot(res_standard$vectors[, 2]*sign(res_standard$vectors[sp_card+1, 2]), main = "Second Eigenvector", xlab = "Index", ylab = "", type = "h")
  lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
  plot(res_standard$vectors[, 3]*sign(res_standard$vectors[2*sp_card+1, 3]), main = "Third Eigenvector", xlab = "Index", ylab = "", type = "h")
  lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")
} else {
  par(mfcol = c(3, 3))
  plot(res_sparse$vectors[, 1]*sign(res_sparse$vectors[1, 1]), main = "First Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
  lines(V[, 1]*sign(V[1, 1]), col = "red")
  plot(res_sparse$vectors[, 2]*sign(res_sparse$vectors[sp_card+1, 2]), main = "Second Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
  lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
  plot(res_sparse$vectors[, 3]*sign(res_sparse$vectors[2*sp_card+1, 3]), main = "Third Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
  lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")

  plot(res_sparseCov$vectors[, 1]*sign(res_sparseCov$vectors[1, 1]), main = "First Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
  lines(V[, 1]*sign(V[1, 1]), col = "red")
  plot(res_sparseCov$vectors[, 2]*sign(res_sparseCov$vectors[sp_card+1, 2]), main = "Second Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
  lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
  plot(res_sparseCov$vectors[, 3]*sign(res_sparseCov$vectors[2*sp_card+1, 3]), main = "Third Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
  lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")

  plot(res_standard$vectors[, 1]*sign(res_standard$vectors[1, 1]), main = "First Eigenvector", xlab = "Index", ylab = "", type = "h")
  lines(V[, 1]*sign(V[1, 1]), col = "red")
  plot(res_standard$vectors[, 2]*sign(res_standard$vectors[sp_card+1, 2]), main = "Second Eigenvector", xlab = "Index", ylab = "", type = "h")
  lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
  plot(res_standard$vectors[, 3]*sign(res_standard$vectors[2*sp_card+1, 3]), main = "Third Eigenvector", xlab = "Index", ylab = "", type = "h")
  lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")
}

