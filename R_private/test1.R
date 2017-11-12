#
# This function is for us to test the package and so
#


source("spEigen.R")

#--------------------#
# Libraries required #
library(mvtnorm) # rmvnorm function for data generation
library(gmodels) # fast.svd function


#------------#
# Parameters #
m = 500 # dimension
n = 100 # number of samples
q = 3 # number of sparse eigenvectors to be estimated
SpCard = 0.2*m # cardinality of the sparse eigenvectors
d = rep(1, q) # IMRP parameter


#-----------------#
# True Covariance #
# Sparse eigenvectors
V = matrix(rnorm(m^2), ncol = m)
tmp = matrix(0, m, q)

for (i in 1:max(q,2)) {
  ind1 = (i-1)*SpCard + 1
  ind2 = i*SpCard
  tmp[ind1:ind2,i] = 1/sqrt(SpCard)
  V[,i] = tmp[,i]
}

V = qr.Q(qr(V)) # orthogonalization, but keep the first eigenvectors to be same as V

# Eigenvalues
vl = rep(1,m)
for (i in 1:q) {
  vl[i] = 100*(q+1-i)
}

# Covariance matrix
R = V %*% diag(vl) %*% t(V)


#-------------#
# Data Matrix #
C = rmvnorm(n = n, mean = rep(0,m), sigma = R) # random data with underlying sparse structure
C = C - matrix(rep(colMeans(C), n), nrow = n, byrow = T) # center the data


#-------------------------------#
# Sparse Eigenvector Extraction #
rho = 0.6
Res = spEigen(C, d, rho, q)


#-------#
# Plots #
Rec = abs(diag(t(Res$sp.vectors) %*% V[,1:q])) # recovery
print(Rec)

par(mfrow = c(3,1))
plot(Res$sp.vectors[,1]*sign(Res$sp.vectors[1,1]), main = "First Sparse Eigenvector", xlab="Index", ylab = "", type = "h")
lines(V[,1]*sign(V[1,1]), col = "red")
plot(Res$sp.vectors[,2]*sign(Res$sp.vectors[SpCard+1,2]), main = "Second Sparse Eigenvector", xlab="Index", ylab = "", type = "h")
lines(V[,2]*sign(V[SpCard+1,2]), col = "red")
plot(Res$sp.vectors[,3]*sign(Res$sp.vectors[2*SpCard+1,3]), main = "Third Sparse Eigenvector", xlab="Index", ylab = "", type = "h")
lines(V[,3]*sign(V[2*SpCard+1,3]), col = "red")

par(mfrow = c(3,1))
plot(Res$vectors[,1]*sign(Res$vectors[1,1]), main = "First Eigenvector", xlab="Index", ylab = "", type = "h")
lines(V[,1]*sign(V[1,1]), col = "red")
plot(Res$vectors[,2]*sign(Res$vectors[SpCard+1,2]), main = "Second Eigenvector", xlab="Index", ylab = "", type = "h")
lines(V[,2]*sign(V[SpCard+1,2]), col = "red")
plot(Res$vectors[,3]*sign(Res$vectors[2*SpCard+1,3]), main = "Third Eigenvector", xlab="Index", ylab = "", type = "h")
lines(V[,3]*sign(V[2*SpCard+1,3]), col = "red")


