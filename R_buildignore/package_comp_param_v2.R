# Libraries required #
library(sparseEigen)
library(mvtnorm) # rmvnorm function for data generation
library(elasticnet)
library(rrcovHD)
library(scales)

set.seed(42)

# Parameters
m <- 500 # dimension
n <- 200 # number of samples
q <- 3 # number of sparse eigenvectors to be estimated
sp_card <- 0.1*m # cardinality of the sparse eigenvectors
rho <- c(0.1, 0.5, 0.9) # sparsity parameters for spEigen()
rho_spca <- c(0.1, 0.5, 0.9) # sparsity parameters for spca()
lam <- c(0.1, 0.5, 0.9) # sparsity parameters for SpcaGrid()

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

# Data Matrix
X <- rmvnorm(n = n, mean = rep(0, m), sigma = R) # random data with underlying sparse structure

# Sparse Eigenvector Extraction
results_spEigen <- array(0, dim=c(3, m, 3))
results_spca <- array(0, dim=c(3, m, 3))
results_spcagrid <- array(0, dim=c(3, m, 3))

angleRec_spEigen <- matrix(0, 3, 3)
angleRec_spca <- matrix(0, 3, 3)
angleRec_spcagrid <- matrix(0, 3, 3)

patternRec_spEigen <- matrix(0, 3, 3)
patternRec_spca <- matrix(0, 3, 3)
patternRec_spcagrid <- matrix(0, 3, 3)

for (i in 1:3) {
  # spEigen()
  res_spEigen <- spEigen(X, q, rho[i], data=TRUE)
  results_spEigen[i, ,] <- matrix(c(res_spEigen$vectors[, 1]*sign(res_spEigen$vectors[1, 1]),
                                    res_spEigen$vectors[, 2]*sign(res_spEigen$vectors[sp_card+1, 2]),
                                    res_spEigen$vectors[, 3]*sign(res_spEigen$vectors[2*sp_card+1, 3])), ncol=3)

  angleRec_spEigen[i, ] <- abs(diag(Conj(t(results_spEigen[i, ,])) %*% V[, 1:q])) # recovery
  patternRec_spEigen[i, ] <- 1 - 1/m * colSums((abs(results_spEigen[i, ,]) > 1e-6) -  (abs(V[, 1:q]) > 1e-6))

  # spca()
  res_spca <- spca(cov(X), K=q, type="Gram", sparse="penalty", trace=FALSE, para=rho_spca[i]*rep(1, 3))
  results_spca[i, ,] <- matrix(c(res_spca$loadings[, 1]*sign(res_spca$loadings[1, 1]),
                                 res_spca$loadings[, 2]*sign(res_spca$loadings[sp_card+1, 2]),
                                 res_spca$loadings[, 3]*sign(res_spca$loadings[2*sp_card+1, 3])), ncol=3)

  angleRec_spca[i, ] <- abs(diag(Conj(t(results_spca[i, ,])) %*% V[, 1:q]))
  patternRec_spca[i, ] <- 1 - 1/m * colSums((abs(results_spca[i, ,]) > 1e-6) -  (abs(V[, 1:q]) > 1e-6))

  # SPcaGrid()
  res_spcagrid <- SPcaGrid(X, lambda=lam[i], k=q)
  results_spcagrid[i, ,] <- matrix(c(getLoadings(res_spcagrid)[, 1]*sign(getLoadings(res_spcagrid)[1, 1]),
                                     getLoadings(res_spcagrid)[, 2]*sign(getLoadings(res_spcagrid)[sp_card+1, 2]),
                                     getLoadings(res_spcagrid)[, 3]*sign(getLoadings(res_spcagrid)[2*sp_card+1, 3])), ncol=3)

  angleRec_spcagrid[i, ] <- abs(diag(Conj(t(results_spcagrid[i, ,])) %*% V[, 1:q]))
  patternRec_spcagrid[i, ] <- 1 - 1/m * colSums((abs(results_spcagrid[i, ,]) > 1e-6) -  (abs(V[, 1:q]) > 1e-6))

}


# Plots
#load("running_time.RData")
cairo_ps(filename = "recovery.ps")
# png(file="recovery.png", width = 20, height = 14, units = "cm", res = 900)
par(mfcol = c(3, 3), mai = c(0.6, 0.3, 0.3, 0.3))
matplot(seq(1, m), results_spEigen[1, ,] , xlab = "", ylab = "", ylim = c(-0.1, 0.25), lty = 'solid',
        main = paste('spEigen: rho =', rho[1]), type = "h",
        col = alpha(c('orangered', 'blue', 'green'), 0.5))
grid()
text(x = 190, y = 0.22, labels = 'Recovery:')
text(x = 300, y = 0.22, labels = 'Angle')
text(x = 300, y = 0.15, labels = round(angleRec_spEigen[1, 1], 4))
text(x = 300, y = 0.10, labels = round(angleRec_spEigen[1, 2], 4))
text(x = 300, y = 0.05, labels = round(angleRec_spEigen[1, 3], 4))
text(x = 400, y = 0.22, labels = 'Pattern')
text(x = 400, y = 0.15, labels = round(patternRec_spEigen[1, 1], 4))
text(x = 400, y = 0.10, labels = round(patternRec_spEigen[1, 2], 4))
text(x = 400, y = 0.05, labels = round(patternRec_spEigen[1, 3], 4))
lines(V[, 1]*sign(V[1, 1]), col = "red")
lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")

matplot(seq(1, m), results_spEigen[2, ,] , xlab = "", ylab = "", ylim = c(-0.1, 0.25), lty = 'solid',
        main = paste('spEigen: rho =', rho[2]), type = "h",
        col = alpha(c('orangered', 'blue', 'green'), 0.5))
grid()
text(x = 190, y = 0.22, labels = 'Recovery:')
text(x = 300, y = 0.22, labels = 'Angle')
text(x = 300, y = 0.15, labels = round(angleRec_spEigen[2, 1], 4))
text(x = 300, y = 0.10, labels = round(angleRec_spEigen[2, 2], 4))
text(x = 300, y = 0.05, labels = round(angleRec_spEigen[2, 3], 4))
text(x = 400, y = 0.22, labels = 'Pattern')
text(x = 400, y = 0.15, labels = round(patternRec_spEigen[2, 1], 4))
text(x = 400, y = 0.10, labels = round(patternRec_spEigen[2, 2], 4))
text(x = 400, y = 0.05, labels = round(patternRec_spEigen[2, 3], 4))
lines(V[, 1]*sign(V[1, 1]), col = "red")
lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")

matplot(seq(1, m), results_spEigen[3, ,] , xlab = "Index", ylab = "", ylim = c(-0.1, 0.25), lty = 'solid',
        main = paste('spEigen: rho =', rho[3]), type = "h",
        col = alpha(c('orangered', 'blue', 'green'), 0.5))
grid()
text(x = 190, y = 0.22, labels = 'Recovery:')
text(x = 300, y = 0.22, labels = 'Angle')
text(x = 300, y = 0.15, labels = round(angleRec_spEigen[3, 1], 4))
text(x = 300, y = 0.10, labels = round(angleRec_spEigen[3, 2], 4))
text(x = 300, y = 0.05, labels = round(angleRec_spEigen[3, 3], 4))
text(x = 400, y = 0.22, labels = 'Pattern')
text(x = 400, y = 0.15, labels = round(patternRec_spEigen[3, 1], 4))
text(x = 400, y = 0.10, labels = round(patternRec_spEigen[3, 2], 4))
text(x = 400, y = 0.05, labels = round(patternRec_spEigen[3, 3], 4))
lines(V[, 1]*sign(V[1, 1]), col = "red")
lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")

matplot(seq(1, m), results_spca[1, ,] , xlab = "", ylab = "", ylim = c(-0.1, 0.25), lty = 'solid',
        main = paste0('spca: para = [', rho_spca[1],', ', rho_spca[1],', ', rho_spca[1],']'), type = "h",
        col = alpha(c('orangered', 'blue', 'green'), 0.5))
grid()
text(x = 190, y = 0.22, labels = 'Recovery:')
text(x = 300, y = 0.22, labels = 'Angle')
text(x = 300, y = 0.15, labels = round(angleRec_spca[1, 1], 4))
text(x = 300, y = 0.10, labels = round(angleRec_spca[1, 2], 4))
text(x = 300, y = 0.05, labels = round(angleRec_spca[1, 3], 4))
text(x = 400, y = 0.22, labels = 'Pattern')
text(x = 400, y = 0.15, labels = round(patternRec_spca[1, 1], 4))
text(x = 400, y = 0.10, labels = round(patternRec_spca[1, 2], 4))
text(x = 400, y = 0.05, labels = round(patternRec_spca[1, 3], 4))
lines(V[, 1]*sign(V[1, 1]), col = "red")
lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")

matplot(seq(1, m), results_spca[2, ,] , xlab = "", ylab = "", ylim = c(-0.1, 0.25), lty = 'solid',
        main = paste0('spca: para = [', rho_spca[2],', ', rho_spca[2],', ', rho_spca[2],']'), type = "h",
        col = alpha(c('orangered', 'blue', 'green'), 0.5))
grid()
text(x = 190, y = 0.22, labels = 'Recovery:')
text(x = 300, y = 0.22, labels = 'Angle')
text(x = 300, y = 0.15, labels = round(angleRec_spca[2, 1], 4))
text(x = 300, y = 0.10, labels = round(angleRec_spca[2, 2], 4))
text(x = 300, y = 0.05, labels = round(angleRec_spca[2, 3], 4))
text(x = 400, y = 0.22, labels = 'Pattern')
text(x = 400, y = 0.15, labels = round(patternRec_spca[2, 1], 4))
text(x = 400, y = 0.10, labels = round(patternRec_spca[2, 2], 4))
text(x = 400, y = 0.05, labels = round(patternRec_spca[2, 3], 4))
lines(V[, 1]*sign(V[1, 1]), col = "red")
lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")

matplot(seq(1, m), results_spca[3, ,] , xlab = "Index", ylab = "", ylim = c(-0.1, 0.25), lty = 'solid',
        main = paste0('spca: para = [', rho_spca[3],', ', rho_spca[3],', ', rho_spca[3],']'), type = "h",
        col = alpha(c('orangered', 'blue', 'green'), 0.5))
grid()
text(x = 190, y = 0.22, labels = 'Recovery:')
text(x = 300, y = 0.22, labels = 'Angle')
text(x = 300, y = 0.15, labels = round(angleRec_spca[3, 1], 4))
text(x = 300, y = 0.10, labels = round(angleRec_spca[3, 2], 4))
text(x = 300, y = 0.05, labels = round(angleRec_spca[3, 3], 4))
text(x = 400, y = 0.22, labels = 'Pattern')
text(x = 400, y = 0.15, labels = round(patternRec_spca[3, 1], 4))
text(x = 400, y = 0.10, labels = round(patternRec_spca[3, 2], 4))
text(x = 400, y = 0.05, labels = round(patternRec_spca[3, 3], 4))
lines(V[, 1]*sign(V[1, 1]), col = "red")
lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")

matplot(seq(1, m), results_spca[1, ,] , xlab = "", ylab = "", ylim = c(-0.1, 0.25), lty = 'solid',
        main = paste0('SPcaGrid: lambda = ', lam[1]), type = "h",
        col = alpha(c('orangered', 'blue', 'green'), 0.5))
grid()
text(x = 190, y = 0.22, labels = 'Recovery:')
text(x = 300, y = 0.22, labels = 'Angle')
text(x = 300, y = 0.15, labels = round(angleRec_spcagrid[1, 1], 4))
text(x = 300, y = 0.10, labels = round(angleRec_spcagrid[1, 2], 4))
text(x = 300, y = 0.05, labels = round(angleRec_spcagrid[1, 3], 4))
text(x = 400, y = 0.22, labels = 'Pattern')
text(x = 400, y = 0.15, labels = round(patternRec_spcagrid[1, 1], 4))
text(x = 400, y = 0.10, labels = round(patternRec_spcagrid[1, 2], 4))
text(x = 400, y = 0.05, labels = round(patternRec_spcagrid[1, 3], 4))
lines(V[, 1]*sign(V[1, 1]), col = "red")
lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")

matplot(seq(1, m), results_spcagrid[2, ,] , xlab = "", ylab = "", ylim = c(-0.1, 0.25), lty = 'solid',
        main = paste0('SPcaGrid: lambda = ', lam[2]), type = "h",
        col = alpha(c('orangered', 'blue', 'green'), 0.5))
grid()
text(x = 190, y = 0.22, labels = 'Recovery:')
text(x = 300, y = 0.22, labels = 'Angle')
text(x = 300, y = 0.15, labels = round(angleRec_spcagrid[2, 1], 4))
text(x = 300, y = 0.10, labels = round(angleRec_spcagrid[2, 2], 4))
text(x = 300, y = 0.05, labels = round(angleRec_spcagrid[2, 3], 4))
text(x = 400, y = 0.22, labels = 'Pattern')
text(x = 400, y = 0.15, labels = round(patternRec_spcagrid[2, 1], 4))
text(x = 400, y = 0.10, labels = round(patternRec_spcagrid[2, 2], 4))
text(x = 400, y = 0.05, labels = round(patternRec_spcagrid[2, 3], 4))
lines(V[, 1]*sign(V[1, 1]), col = "red")
lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")

matplot(seq(1, m), results_spcagrid[3, ,] , xlab = "Index", ylab = "", ylim = c(-0.1, 0.25), lty = 'solid',
        main = paste0('SPcaGrid: lambda = ', lam[3]), type = "h",
        col = alpha(c('orangered', 'blue', 'green'), 0.5))
grid()
text(x = 190, y = 0.22, labels = 'Recovery:')
text(x = 300, y = 0.22, labels = 'Angle')
text(x = 300, y = 0.15, labels = round(angleRec_spcagrid[3, 1], 4))
text(x = 300, y = 0.10, labels = round(angleRec_spcagrid[3, 2], 4))
text(x = 300, y = 0.05, labels = round(angleRec_spcagrid[3, 3], 4))
text(x = 400, y = 0.22, labels = 'Pattern')
text(x = 400, y = 0.15, labels = round(patternRec_spcagrid[3, 1], 4))
text(x = 400, y = 0.10, labels = round(patternRec_spcagrid[3, 2], 4))
text(x = 400, y = 0.05, labels = round(patternRec_spcagrid[3, 3], 4))
lines(V[, 1]*sign(V[1, 1]), col = "red")
lines(V[, 2]*sign(V[sp_card+1, 2]), col = "red")
lines(V[, 3]*sign(V[2*sp_card+1, 3]), col = "red")

dev.off()
