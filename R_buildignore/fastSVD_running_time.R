# Libraries required
library(gmodels) # fast.svd function

# Parameters
m <- seq(100, 1000, 100)  # dimension
n <- c(0.05, 0.20, 0.50, 1.0)  # data as a percentage of the dimension
mc <- 200  # monte carlo


##### SVD Vs FAST.SVD (Data)#####
time_svd <- matrix(rep(0, length(m)*length(n)), nrow = length(m))
time_fast_svd <- matrix(rep(0, length(m)*length(n)), nrow = length(m))


for (ii in 1:length(m)) {
  print(m[ii])
  for (jj in 1:length(n)) {
    print(n[jj])
    for (kk in 1:mc) {
      X <- matrix( rnorm(m[ii]*ceiling(n[jj]*m[ii]), mean=0, sd=2), ceiling(n[jj]*m[ii]), m[ii])

      ptm <- proc.time()
      res1 <- svd(X)
      time_svd[ii, jj] <- time_svd[ii, jj] + ( (proc.time() - ptm)[3] ) / mc

      ptm <- proc.time()
      res2 <- fast.svd(X)
      time_fast_svd[ii, jj] <- time_fast_svd[ii, jj] + ( (proc.time() - ptm)[3] ) / mc
    }
  }
}

results <- matrix(cbind(time_svd, time_fast_svd), ncol = 8)

par(mfrow=c(1,1))
matplot(m, results, pch=1, col = 1:8, type = 'b', xlab = "Dimension", ylab = "Time", log='y')
legend("topleft", legend = c('svd 5%', 'svd 20%', 'svd 50%', 'svd 100%', 'fast.svd 5%', 'fast.svd 20%', 'fast.svd 50%', 'fast.svd 100%'), col=1:8, pch=1)
grid()




##### EVD Vs SVD (Covariance)#####
time_svdCov <- rep(0, length(m))
time_evdCov <- rep(0, length(m))


for (ii in 1:length(m)) {
  print(m[ii])
  for (kk in 1:mc) {
    # Consider samples = 1.5 times dimension (full rank)
    X <- matrix( rnorm(m[ii]*ceiling(1.5*m[ii]), mean=0, sd=2), ceiling(n[jj]*m[ii]), m[ii])
    S = cov(X)

    ptm <- proc.time()
    res1 <- svd(S)
    time_svdCov[ii] <- time_svdCov[ii] + ( (proc.time() - ptm)[3] ) / mc

    ptm <- proc.time()
    res2 <- eigen(S)
    time_evdCov[ii] <- time_evdCov[ii] + ( (proc.time() - ptm)[3] ) / mc
  }
}


resultsCov <- matrix(cbind(time_svdCov, time_evdCov), ncol = 2)

matplot(m, resultsCov, pch=1, col = 1:2, type = 'b', xlab = "Dimension", ylab = "Time")
legend("topleft", legend = c('svd', 'evd'), col=1:2, pch=1)
grid()






