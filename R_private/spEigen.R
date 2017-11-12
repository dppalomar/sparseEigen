spEigen = function(C, d, rho_nrm, q, V = NULL){

# INPUT
#   C :             n-by-m data matrix (n samples, m variables).
#   d :             1-by-q vector with weights.
#   rho_nrm :       Sparsity weight factor. Values from 0 to 1.
#   q :             Number of estimated eigenvectors.
#   V (optional):   m-by-q initial point matrix. If not provided the eigenvectors
#                   of the sample covariance matrix are used.
#
# OUTPUT
#   sp.vectors :    m-by-q matrix, columns corresponding to leading
#                   sparse eigenvectors.
#   vectors :       m-by-q matrix, columns corresponding to leading
#                   eigenvectors.
#   values :        q-by-1 vector corresponding to the leading eigenvalues.
#
# INFO
#   Reference:      K. Benidis, Y. Sun, P. Babu, D.P. Palomar "Orthogonal Sparse
#                   PCA and Covariance Estimation via Procrustes Reformulation"
#                   IEEE Transactions on Signal Processing, vol 64, Dec. 2016.
#
#   Algorithm:      This algorithm corresponds to the accelerated IMRP algorithm
#                   of the referenced paper.
#
#   Link :          http://www.danielppalomar.com/publications.html


# Initialize
m = ncol(C)
n = nrow(C)
k = 0
maxIter = 1000

# Preallocation
V_tld = matrix(0, m, q)
H = matrix(0, m, q)
F_v = matrix(0, maxIter, 1)
g = matrix(0, m, q)

# Rho
svd_c = fast.svd(C)
Sc2 = svd_c$d^2
rho = rho_nrm*max(colSums(C^2))*(Sc2[1:q]/Sc2[1])*d

# Initial Point
if (is.null(V)) {
  V = svd_c$v[,1:q]
}

# Decreasing epsilon, p
K = 10
p1 = 1 # first value of p
pT = 7 # last value of p
gamma = (pT/p1)^(1/K)
pp = p1*gamma^(0:K)
pp = 10^(-pp)

tol = pp*1e-2 # tolerance for convergence
Eps = pp # epsilon

for (ee in 1:(K+1)) {
  p = pp[ee]
  epsi = Eps[ee]
  c1 = log(1 + 1/p)
  c2 = 2*(p + epsi)*c1
  w0 = (1/(epsi*c2))*rep(1,m*q)
  flg = 1

  while (1) {
    k = k + 1

    #-------------------------------------#
    # First iteration of the acceleration #

    # weights
    w = w0
    ind = abs(c(V)) > epsi
    w[ind] = (0.5/c1)/(V[ind]^2 + p*abs(V[ind]))

    # MM
    for (i in 1:q) {
      w_tmp = w[((i-1)*m+1):(i*m)]
      V_tld[,i] = V[,i]*d[i]
      H[,i] = (w_tmp - max(w_tmp)*rep(1,m))*V[,i]*rho[i]
    }

    G = svd_c$v %*% ((t(svd_c$v) %*% V_tld)*matrix(rep(Sc2, q), ncol = q))

    # update
    s1 = fast.svd(G - H)
    V1 = s1$u %*% t(s1$v)

    #--------------------------------------#
    # Second iteration of the acceleration #

    # weights
    w = w0
    ind = abs(c(V1)) > epsi
    w[ind] = (0.5/c1)/(V1[ind]^2 + p*abs(V1[ind]))

    # MM
    for (i in 1:q) {
      w_tmp = w[((i-1)*m+1):(i*m)]
      V_tld[,i] = V1[,i]*d[i]
      H[,i] = (w_tmp - max(w_tmp)*rep(1,m))*V1[,i]*rho[i]
    }

    G = svd_c$v %*% ((t(svd_c$v) %*% V_tld)*matrix(rep(Sc2, q), ncol = q))

    # update
    s2 = fast.svd(G - H)
    V2 = s2$u %*% t(s2$v)

    #--------------#
    # Acceleration #

    R = V1 - V
    U = V2 - V1 - R
    a = min(-norm(R, type = "F")/norm(U, type = "F"), -1)

    while (1) { # backtracking loop
      V0 = V - 2*a*R + a^2*U

      # Projection
      s3 = fast.svd(V0)
      V0 = s3$u %*% t(s3$v)

      g[abs(V0) <= epsi] = V0[abs(V0) <= epsi]^2/(epsi*c2)
      g[abs(V0) > epsi] = log((p + abs(V0[abs(V0) > epsi]))/(p + epsi))/c1 + epsi/c2

      F_v[k] = colSums((C %*% V0)^2) %*% d - colSums(g) %*% rho

      if (flg == 0 && F_v[k]*(1 + sign(F_v[k])*1e-9) <= F_v[max(k-1,1)]) {
        a = (a-1)/2
      }
      else {
        V = V0
        break
      }
    }

    # Stopping criterion
    if (flg == 0) {
      rel_change = abs(F_v[k] - F_v[k-1])/max(1, abs(F_v[k-1])) # relative change in objective
      if (rel_change <= tol[ee] || k >= maxIter) {
        F_v = F_v[1:k]
        break
      }
    }
    flg = 0
  }
}

V[abs(V) < 1e-10] = 0; # threshold
nrm = 1/sqrt(colSums(V^2))
V = matrix(rep(nrm, m), ncol = q)*V

return(list(sp.vectors = V, vectors = svd_c$v[,1:q], values = Sc2[1:q]/n))

}

