source("spEigen.R")

#--------------------#
# Libraries required #
library(mvtnorm) # rmvnorm function for data generation

#------------#
# Parameters #
m <- 500 # dimension
n <- 200 # number of samples
q <- 3 # number of sparse eigenvectors to be estimated
sp_card <- 0.1*m # cardinality of the sparse eigenvectors
rho <- 0.5

#-----------------#
# True Covariance #
# Sparse eigenvectors
V <- matrix(rnorm(m^2)*exp(1i*runif(m^2,0,2*pi)), ncol = m)
tmp <- matrix(0, m, q)

for (i in 1:max(q, 2)) {
  ind1 <- (i-1)*sp_card + 1
  ind2 <- i*sp_card
  tmp[ind1:ind2, i] = 1/sqrt(sp_card)
  V[, i] <- tmp[, i]*exp(1i*runif(m,0,2*pi))
}
for (i in (q+1):m) {
  tmp <- V[,i]
  for (j in 1:(i-1)) {
    tmp <- tmp - (V[,i] %*% Conj(V[,j])) %*% V[,j]
  }
  V[,i] <- tmp/sqrt(sum(abs(tmp)^2))
}

# Eigenvalues
vl <- rep(1, m)
for (i in 1:q) {
  vl[i] <- 100*(q + 1 - i)
}

# Covariance matrix
R <- V %*% diag(vl) %*% Conj(t(V))

#-------------#
# Data Matrix #
X <- rmvnorm(n = n, mean = rep(0, m), sigma = R) # random data with underlying sparse structure

#-------------------------------#
# Sparse Eigenvector Extraction #
data <- FALSE

if (data) {
  res_sparse <- spEigen(X, q, rho, data = TRUE)
} else {
  X <- scale(X, center = TRUE, scale = FALSE)
  S <- 1/(n-1) * Conj(t(X)) %*% X
  res_sparse <- spEigen(S, q, rho, data = TRUE)
}

#-------#
# Plots #
par(mfcol = c(3, 1))
plot(abs(res_sparse$vectors[, 1]), main = "First Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(abs(V[, 1]), col = "red")
plot(abs(res_sparse$vectors[, 2]), main = "Second Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(abs(V[, 2]), col = "red")
plot(abs(res_sparse$vectors[, 3]), main = "Third Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(abs(V[, 3]), col = "red")

