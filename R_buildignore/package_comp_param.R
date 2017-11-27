source("spEigen.R")

#--------------------#
# Libraries required #
library(mvtnorm) # rmvnorm function for data generation
library(elasticnet)
library(scales)

set.seed(42)

#------------#
# Parameters #
m <- 500 # dimension
n <- 200 # number of samples
q <- 3 # number of sparse eigenvectors to be estimated
sp_card <- 0.1*m # cardinality of the sparse eigenvectors
rho <- c(0.1, 0.5, 0.9)
rho_spca <- c(0.1, 0.5, 0.9)

#-----------------#
# True Covariance #
# generate non-overlapping sparse eigenvectors
V <- matrix(rep(0, m*q), ncol = q)
V[cbind(seq(1, q*sp_card), rep(1:q, each = sp_card))] <- 1/sqrt(sp_card)
V <- cbind(V, matrix(rnorm(m*(m-q)), m, m-q))

V <- qr.Q(qr(V)) # orthogonalization, but keep the first eigenvectors to be same as V

# generate eigenvalues
lmd <- c(100*seq(from = q, to = 1), rep(1, m-q))

# generate covariance matrix from sparse eigenvectors and eigenvalues
R <- V %*% diag(lmd) %*% Conj(t(V))

#-------------#
# Data Matrix #
X <- rmvnorm(n = n, mean = rep(0, m), sigma = R) # random data with underlying sparse structure

#-------------------------------#
# Sparse Eigenvector Extraction #
results_spEigen <- array(0, dim=c(3, m, 3))
results_spca <- array(0, dim=c(3, m, 3))

for (i in 1:3) {
  res_spEigen <- spEigen(X, q, rho[i], data=TRUE)
  results_spEigen[i, ,] <- matrix(c(res_spEigen$vectors[, 1]*sign(res_spEigen$vectors[1, 1]),
                                 res_spEigen$vectors[, 2]*sign(res_spEigen$vectors[sp_card+1, 2]),
                                 res_spEigen$vectors[, 3]*sign(res_spEigen$vectors[2*sp_card+1, 3])), ncol=3)

  res_spca <- spca(cov(X), K=q, type="Gram", sparse="penalty", trace=FALSE, para=rho_spca[i]*rep(1, 3))
  results_spca[i, ,] <- matrix(c(res_spca$loadings[, 1]*sign(res_spca$loadings[1, 1]),
                                 res_spca$loadings[, 2]*sign(res_spca$loadings[sp_card+1, 2]),
                                 res_spca$loadings[, 3]*sign(res_spca$loadings[2*sp_card+1, 3])), ncol=3)
}


#----------#
# Recovery #
recovery_spEigen <- matrix(0, 3, 3)
recovery_spca <- matrix(0, 3, 3)

for (i in 1:3) {
  recovery_spEigen[i, ] <- abs(diag(Conj(t(res_spEigen$vectors)) %*% V[, 1:q])) # recovery
  recovery_spca[i, ] <- abs(diag(Conj(t(res_spca$loadings)) %*% V[, 1:q])) # recovery
}


#-------#
# Plots #
par(mfcol = c(3, 2))
matplot(seq(1, m), results_spEigen[1, ,] , xlab = "", ylab = "", ylim = c(-0.1, 0.2),
        main = paste('spEigen: rho =', rho[1]), type = "h",
        col = alpha(c('orangered', 'blue', 'green'), 0.5))
text(x = 300, y = 0.15, labels = round(recovery_spEigen[1, 1], 4))
text(x = 300, y = 0.10, labels = round(recovery_spEigen[1, 2], 4))
text(x = 300, y = 0.05, labels = round(recovery_spEigen[1, 3], 4))
# legend("topright", legend = c('1st', '2nd', '3rd'), col=c('orangered', 'blue', 'green'), pch=1)
lines(V[, 1]*sign(V[1, 1]), col = "red")
lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")
grid()

matplot(seq(1, m), results_spEigen[2, ,] , xlab = "", ylab = "", ylim = c(-0.1, 0.2),
        main = paste('spEigen: rho =', rho[2]), type = "h",
        col = alpha(c('orangered', 'blue', 'green'), 0.5))
text(x = 300, y = 0.15, labels = round(recovery_spEigen[2, 1], 4))
text(x = 300, y = 0.10, labels = round(recovery_spEigen[2, 2], 4))
text(x = 300, y = 0.05, labels = round(recovery_spEigen[2, 3], 4))
lines(V[, 1]*sign(V[1, 1]), col = "red")
lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")
grid()

matplot(seq(1, m), results_spEigen[3, ,] , xlab = "Index", ylab = "", ylim=c(-0.1, 0.2),
        main = paste('spEigen: rho =', rho[3]), type = "h",
        col = alpha(c('orangered', 'blue', 'green'), 0.5))
text(x = 300, y = 0.15, labels = round(recovery_spEigen[3, 1], 4))
text(x = 300, y = 0.10, labels = round(recovery_spEigen[3, 2], 4))
text(x = 300, y = 0.05, labels = round(recovery_spEigen[3, 3], 4))
lines(V[, 1]*sign(V[1, 1]), col = "red")
lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")
grid()

matplot(seq(1, m), results_spca[1, ,] , xlab = "", ylab = "", ylim=c(-0.1, 0.2),
        main = paste0('spca: para = [', rho_spca[1],', ', rho_spca[1],', ', rho_spca[1],']'), type = "h",
        col = alpha(c('orangered', 'blue', 'green'), 0.5))
text(x = 300, y = 0.15, labels = round(recovery_spca[1, 1], 4))
text(x = 300, y = 0.10, labels = round(recovery_spca[1, 2], 4))
text(x = 300, y = 0.05, labels = round(recovery_spca[1, 3], 4))
lines(V[, 1]*sign(V[1, 1]), col = "red")
lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")
grid()

matplot(seq(1, m), results_spca[2, ,] , xlab = "", ylab = "", ylim=c(-0.1, 0.2),
        main = paste0('spca: para = [', rho_spca[2],', ', rho_spca[2],', ', rho_spca[2],']'), type = "h",
        col = alpha(c('orangered', 'blue', 'green'), 0.5))
text(x = 300, y = 0.15, labels = round(recovery_spca[2, 1], 4))
text(x = 300, y = 0.10, labels = round(recovery_spca[2, 2], 4))
text(x = 300, y = 0.05, labels = round(recovery_spca[2, 3], 4))
lines(V[, 1]*sign(V[1, 1]), col = "red")
lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")
grid()

matplot(seq(1, m), results_spca[3, ,] , xlab = "Index", ylab = "", ylim=c(-0.1, 0.2),
        main = paste0('spca: para = [', rho_spca[3],', ', rho_spca[3],', ', rho_spca[3],']'), type = "h",
        col = alpha(c('orangered', 'blue', 'green'), 0.5))
text(x = 300, y = 0.15, labels = round(recovery_spca[3, 1], 4))
text(x = 300, y = 0.10, labels = round(recovery_spca[3, 2], 4))
text(x = 300, y = 0.05, labels = round(recovery_spca[3, 3], 4))
lines(V[, 1]*sign(V[1, 1]), col = "red")
lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")
grid()

