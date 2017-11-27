source("spEigen.R")
source("spEigenCov.R")

#--------------------#
# Libraries required #
library(cmvnorm) # rmvnorm function for data generation

#------------#
# Parameters #
m <- 500 # dimension
n <- 600 # number of samples
q <- 3 # number of sparse eigenvectors to be estimated
sp_card <- 0.2*m # cardinality of the sparse eigenvectors
rho <- 0.5

# generate non-overlapping sparse eigenvectors
V <- matrix(rep(0, m*q), ncol = q)
V[cbind(seq(1, q*sp_card), rep(1:q, each = sp_card))] <- 1/sqrt(sp_card) * exp(1i*runif(sp_card*q, 0, 2*pi))
V <- cbind(V, matrix(rnorm(m*(m-q))*exp(1i*runif(m*(m-q),0,2*pi)), m, m-q))

# keep first q eigenvectors the same (already orthogonal) and orthogonalize the rest
for (i in (q+1):m) {
  tmp <- V[,i]
  for (j in 1:(i-1)) {
    tmp <- tmp - (V[,i] %*% Conj(V[,j])) %*% V[,j]
  }
  V[,i] <- tmp/sqrt(sum(abs(tmp)^2))
}

# generate eigenvalues
lmd <- c(100*seq(from = q, to = 1), rep(1, m-q))

# generate covariance matrix from sparse eigenvectors and eigenvalues
R <- V %*% diag(lmd) %*% Conj(t(V))

#-------------#
# Data Matrix #
X <- rcmvnorm(n = n, mean = rep(0, m), sigma = R) # random data with underlying sparse structure
X <- scale(X, center = TRUE, scale = FALSE)

#-------------------------------#
# Sparse Eigenvector Extraction #
data <- FALSE

if (data) {
  S <- 1/(n-1) * Conj(t(X)) %*% X
  res_sparse <- spEigen(X, q, rho, data = TRUE)
  res_sparseCov <- spEigenCov(S, q, rho)
} else {
  S <- 1/(n-1) * Conj(t(X)) %*% X
  res_sparse <- spEigen(S, q, rho)
  res_sparseCov <- spEigenCov(S, q, rho)
}

#-------#
# Plots #
par(mfcol = c(3, 2))
plot(abs(res_sparse$vectors[, 1]), main = "First Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(abs(V[, 1]), col = "red")
plot(abs(res_sparse$vectors[, 2]), main = "Second Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(abs(V[, 2]), col = "red")
plot(abs(res_sparse$vectors[, 3]), main = "Third Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(abs(V[, 3]), col = "red")

plot(abs(res_sparseCov$vectors[, 1]), main = "First Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(abs(V[, 1]), col = "red")
plot(abs(res_sparseCov$vectors[, 2]), main = "Second Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(abs(V[, 2]), col = "red")
plot(abs(res_sparseCov$vectors[, 3]), main = "Third Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(abs(V[, 3]), col = "red")

# show error between estimated and true covariance
norm(Re(S) - Re(R), type = 'F') #for sample covariance matrix
norm(Re(res_sparseCov$cov) - Re(R), type = 'F') #for covariance with sparse eigenvectors

norm(abs(Im(S)) - abs(Im(R)), type = 'F') #for sample covariance matrix
norm(abs(Im(res_sparseCov$cov)) - abs(Im(R)), type = 'F') #for covariance with sparse eigenvectors

norm(abs(S) - abs(R), type = 'F')
norm(abs(res_sparseCov$cov) - abs(R), type = 'F')
