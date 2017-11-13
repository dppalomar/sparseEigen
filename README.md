<!-- README.md is generated from README.Rmd. Please edit that file -->
sparseEigen
===========

The goal of sparseEigen is to estimate sparse eigenvectors of numerical matrices. The orthogonality property of the eigenvectors is retained.

Installation
------------

``` r
# Installation from local file
install.packages(file.choose(), repos = NULL, type="source")

# Or from GitHub
# install.packages("devtools")
devtools::install_github("dppalomar/sparseEigen")

# Get help
library(sparseEigen)
help(package="sparseEigen")
?spEigen
```

Example
-------

This is a simple illustrative example:

``` r
# PARAMETERS 
m <- 500 # dimension
n <- 100 # number of samples
q <- 3 # number of sparse eigenvectors to be estimated
SpCard <- 0.2*m # cardinality of the sparse eigenvectors
rho <- 0.6 # sparsity level

# DATA
# Sparse artificial eigenvectors
V <- matrix(rnorm(m^2), ncol = m)
tmp <- matrix(0, m, q)

for (i in 1:max(q, 2)) {
  ind1 <- (i-1)*SpCard + 1
  ind2 <- i*SpCard
  tmp[ind1:ind2, i] = 1/sqrt(SpCard)
  V[, i] <- tmp[, i]
}

V <- qr.Q(qr(V)) # orthogonalization, but keep the first eigenvectors to be same as V

# Eigenvalues
vl <- rep(1, m)
for (i in 1:q) {
  vl[i] <- 100*(q + 1 - i)
}

# Covariance matrix from sparse eigenvectors and eigenvalues
R <- V %*% diag(vl) %*% t(V)

# Data matrix - Generate data from a zero mean multivariate Gaussian and the constructed covariance
C <- mvtnorm::rmvnorm(n = n, mean = rep(0, m), sigma = R) # random data with underlying sparse structure
C <- C - matrix(rep(colMeans(C), n), nrow = n, byrow = T) # center the data

# ALGORITHM
res <- spEigen(C, q, rho)

# EVALUATION
# Recovery
recovery <- abs(diag(t(res$sp.vectors) %*% V[, 1:q])) # recovery
print(recovery)

# Comparison to normal eigenvectors
par(mfcol = c(3, 2))
plot(res$sp.vectors[, 1]*sign(res$sp.vectors[1, 1]), main = "First Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(V[, 1]*sign(V[1, 1]), col = "red")
plot(res$sp.vectors[, 2]*sign(res$sp.vectors[SpCard+1, 2]), main = "Second Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(V[, 2]*sign(V[SpCard+1, 2]), col = "red")
plot(res$sp.vectors[, 3]*sign(res$sp.vectors[2*SpCard+1, 3]), main = "Third Sparse Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(V[, 3]*sign(V[2*SpCard+1, 3]), col = "red")

plot(res$vectors[, 1]*sign(res$vectors[1, 1]), main = "First Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(V[, 1]*sign(V[1, 1]), col = "red")
plot(res$vectors[, 2]*sign(res$vectors[SpCard+1, 2]), main = "Second Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(V[, 2]*sign(V[SpCard+1, 2]), col = "red")
plot(res$vectors[, 3]*sign(res$vectors[2*SpCard+1, 3]), main = "Third Eigenvector", xlab = "Index", ylab = "", type = "h")
lines(V[, 3]*sign(V[2*SpCard+1, 3]), col = "red")
```
