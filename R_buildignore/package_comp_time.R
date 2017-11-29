
# Libraries required
library(sparseEigen)
library(mvtnorm) # rmvnorm function for data generation
library(elasticnet)

# Constant parameters
q <- 3 # number of sparse eigenvectors to be estimated
rho <- 0.5 # sparsity
mc <- 50 # monte carlo for each dimension


########## Loop for different dimensions ##########
dims <- 100*seq(1, 8)

# Initialize
time_spEigen <- matrix(rep(0, mc*length(dims)), nrow = length(dims))
time_spca <- matrix(rep(0, mc*length(dims)), nrow = length(dims))


for (dd in 1:length(dims)) {
  print(dims[dd])

  # Parameters that change in every loop
  m <- dims[dd] # dimension
  n <- ceiling(m/5) # number of samples
  sp_card <- ceiling(0.1*m) # cardinality of the sparse eigenvectors

  ### True Covariance ###
  # generate non-overlapping sparse eigenvectors
  V <- matrix(rep(0, m*q), ncol = q)
  V[cbind(seq(1, q*sp_card), rep(1:q, each = sp_card))] <- 1/sqrt(sp_card)
  V <- cbind(V, matrix(rnorm(m*(m-q)), m, m-q))

  V <- qr.Q(qr(V)) # orthogonalization, but keep the first eigenvectors to be same as V

  # generate eigenvalues
  lmd <- c(100*seq(from = q, to = 1), rep(1, m-q))

  # generate covariance matrix from sparse eigenvectors and eigenvalues
  R <- V %*% diag(lmd) %*% Conj(t(V))

  # Call the function once to get it into memory
  res_spEigen <- spEigen(R, q, rho)
  res_spca <- spca(R, K=q, type="Gram", sparse="penalty", trace=FALSE, para=c(.4,.4,.4))

  ########## Monte Carlo ##########
  for (ii in 1:mc) {
    # Data Matrix
    X <- rmvnorm(n = n, mean = rep(0, m), sigma = R) # random data with underlying sparse structure

    # 1. spEigen()
    ptm <- proc.time()
    res_spEigen <- spEigen(X, q, rho, data=TRUE)
    time_spEigen[dd, ii] <- (proc.time() - ptm)[3]

    # 2. spca() - without including the time of calculating the covariance
    Xc <- cov(X)
    ptm <- proc.time()
    res_spca <- spca(Xc, K=q, type="Gram", sparse="penalty", trace=FALSE, para=c(.4,.4,.4))
    time_spca[dd, ii] <- (proc.time() - ptm)[3]
  }
}

# average
avg_spEigen <- rowMeans(time_spEigen)
avg_spca <- rowMeans(time_spca)
results <- matrix(c(avg_spca, avg_spEigen), ncol=2)

########## Plots ##########
#load("running_time.RData")
png(file="running_time2.png", width = 18, height = 13, units = "cm", res = 600)
matplot(dims, results, pch=1, col = 1:2, type = 'b',
        xlab = "Dimension", ylab = "Time", log = 'y', yaxt = 'n')
axis(2, at = 10^(c(-1, 0, 1, 2)))
legend("topleft", legend = c('spca()', 'spEigen()'), col=1:2, pch=1)
grid()
par(new = TRUE)
matplot(dims, results, pch=1, col = 1:2, type = 'b',
        xlab = "Dimension", ylab = "Time", log = 'y', yaxt = 'n')
axis(2, at = 10^(c(-1, 0, 1, 2)))
legend("topleft", legend = c('spca()', 'spEigen()'), col=1:2, pch=1)
dev.off()



