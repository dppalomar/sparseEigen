# Test the running time in the following cases:
# 1. Use the data matrix
# 2. Use covariance
# 3. Transform the data matrix to covariance
# 4. Transfrom the covariance to 'data' through Cholesky factorization

source("spEigen_test.R")

# Libraries required 
library(mvtnorm) # rmvnorm function for data generation
library(gmodels) # fast.svd function

# Constant parameters
q <- 3 # number of sparse eigenvectors to be estimated
rho <- 0.6 # sparsity 
mc <- 200 # monte carlo for each dimension


########## Loop for different dimensions ##########
dims <- 100*seq(1, 10)

# Initialize
time_data <- matrix(rep(0, mc*length(dims)), nrow = length(dims))
time_cov <- matrix(rep(0, mc*length(dims)), nrow = length(dims))
time_data_transf <- matrix(rep(0, mc*length(dims)), nrow = length(dims))
time_cov_transf <- matrix(rep(0, mc*length(dims)), nrow = length(dims))


for (dd in 1:length(dims)) {
  print(dims[dd])
  
  # Parameters that change in every loop
  m <- dims[dd] # dimension
  n <- ceiling(m/5) # number of samples
  sp_card <- ceiling(0.2*m) # cardinality of the sparse eigenvectors

  # True Covariance 
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

  # Call the function once to get it into memory
  res_sparse <- spEigen_test(R, q, rho)
  
  ########## Monte Carlo ##########
  for (ii in 1:mc) {
    # Data Matrix 
    X <- rmvnorm(n = n, mean = rep(0, m), sigma = R) # random data with underlying sparse structure
    X <- scale(X, center = TRUE, scale = FALSE) # center the data

    # 1. Data matrix
    ptm <- proc.time()
    res_sparse <- spEigen_test(X, q, rho, cov = FALSE)
    time_data[dd, ii] <- (proc.time() - ptm)[3]    
  
    # 2. Covariance matrix
    Xc <- cov(X) # do not count the covariance computation
    ptm <- proc.time()
    res_sparse <- spEigen_test(Xc, q, rho)
    time_cov[dd, ii] <- (proc.time() - ptm)[3]    

    # 3. Data matrix to covariance
    ptm <- proc.time()
    res_sparse <- spEigen_test(X, q, rho, cov = FALSE, transf = TRUE)
    time_data_transf[dd, ii] <- (proc.time() - ptm)[3]    

    # 4. Covariance to data
    ptm <- proc.time()
    res_sparse <- spEigen_test(Xc, q, rho, transf = TRUE)
    time_cov_transf[dd, ii] <- (proc.time() - ptm)[3]    
  }
}

# average
avg_data <- rowMeans(time_data)
avg_cov <- rowMeans(time_cov)
avg_data_transf <- rowMeans(time_data_transf)
avg_cov_transf <- rowMeans(time_cov_transf)
results <- matrix(c(avg_data, avg_cov, avg_data_transf, avg_cov_transf), ncol = 4)

########## Plots ##########
matplot(dims, results , pch=1, col = 1:4, type = 'b', main = "Average running time", xlab = "Dimension", ylab = "Time", log='y') 
legend("topleft", legend = c('Data', 'Covariance', 'Data transformed', 'Covariance transformed'), col=1:4, pch=1) 
grid()




