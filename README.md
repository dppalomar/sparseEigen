<!-- README.md is generated from README.Rmd. Please edit that file -->
sparseEigen
===========

This package provides two functions to compute sparse eigenvectors (while keeping their orthogonality property) from either the covariance matrix or directly the data matrix.

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

Usage
-----

This is a simple illustrative example:

``` r
library(sparseEigen)

# parameters 
m <- 500  # dimension
n <- 100  # number of samples
q <- 3  # number of sparse eigenvectors to be estimated
sp_card <- 0.2*m  # cardinality of each sparse eigenvector
rho <- 0.6  # sparsity level

# generate non-overlapping sparse eigenvectors
V <- matrix(rnorm(m^2), ncol = m)
tmp <- matrix(0, m, q)
for (i in 1:max(q, 2)) {
  ind1 <- (i - 1)*sp_card + 1
  ind2 <- i*sp_card
  tmp[ind1:ind2, i] = 1/sqrt(sp_card)
  V[, i] <- tmp[, i]
}
V <- qr.Q(qr(V))  # keep first q eigenvectors the same (already orthogonal) and orthogonalize the rest

# generate eigenvalues
lmd <- rep(1, m)
lmd[1:q] <- 100*seq(from = q, to = 1)

# generate covariance matrix from sparse eigenvectors and eigenvalues
R <- V %*% diag(lmd) %*% t(V)

# generate data matrix from a zero-mean multivariate Gaussian distribution with the constructed covariance
X <- MASS::mvrnorm(n, rep(0, m), R)  # random data with underlying sparse structure
X <- scale(X, center = TRUE, scale = FALSE)  # center the data

# computation of sparse eigenvectors
res_standard <- eigen(cov(X))
#res_sparse <- spEigen(cov(X), q, rho)
#res_sparse <- spEigenDataMatrix(X, q, rho)

# show inner product between estimated eigenvectors and originals
inner_product_standard <- abs(diag(t(res_standard$vectors) %*% V[, 1:q]))
#inner_product_sparse <- abs(diag(t(res_sparse$vectors) %*% V[, 1:q]))
cat("Inner product for standard estimated eigenvectors: ", inner_product_standard)
#> Inner product for standard estimated eigenvectors:  0.9903066 0.979756 0.9604419
#cat("Inner product for sparse estimated eigenvectors: ", inner_product_sparse)
```

The following plot shows the sparsity pattern of the eigenvectors: ![](README-unnamed-chunk-4-1.png)