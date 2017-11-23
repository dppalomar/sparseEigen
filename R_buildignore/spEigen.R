spEigen <- function(X, q = 1, rho = 0.5, data = FALSE, d = NA, V = NA, thres = 1e-9) {
  max_iter <- 1000 # maximum MM iterations

  ######## error control  #########
  X <- as.matrix(X)
  m <- ncol(X)
  if (m == 1) stop("Data is univariate!")
  if (q > qr(X)$rank) stop("The number of estimated eigenvectors q should not be larger than rank(X).")
  if (anyNA(X) || anyNA(q) || anyNA(rho)) stop("This function cannot handle NAs.")
  if ((q %% 1) != 0 || q <= 0) stop("The input argument q should be a positive integer.")
  if (rho <= 0) stop("The input argument rho should be positive.")
  #################################

  # Center the data matrix (if data = TRUE)
  if (data)
    X <- scale(X, center = TRUE, scale = FALSE)

  # Input parameter d: vector of weights
  if (is.na(d))
    d <- seq(from = 1, to = 0.5, length.out = q)

  # Sparsity parameter rho
  if (data) {
    svd_x <- svd(X)
    sv2 <- svd_x$d^2
    Vx <- svd_x$v
    rho <- rho * max(colSums(abs(X)^2)) * (sv2[1:q]/sv2[1]) * d
  }
  else {
    eig_x <- eigen(X)
    sv2 <- eig_x$values
    Vx <- eig_x$vectors
    rho <- rho * max(Re(diag(X))) * (sv2[1:q]/sv2[1]) * d
  }

  # Input parameter V: initial point
  if (is.na(V))
    V <- Vx[, 1:q]

  # Preallocation
  F_v <- matrix(0, max_iter, 1)  # objective value
  g <- matrix(0, m, q)

  # Decreasing epsilon, p
  K <- 10  # number of outer iterations
  p1 <- 1  # first value of -log(p)
  pk <- 7  # last value of -log(p)
  gamma <- (pk/p1)^(1/K)
  pp <- p1 * gamma^(0:K)
  pp <- 10^(-pp)

  tol <- pp * 1e-2  # tolerance for convergence
  Eps <- pp  # epsilon


  ######################### MM LOOP #########################
  k <- 0  # MM iteration counter
  for (ee in 1:(K+1)) {  # loop for approximation based on p & epsilon
    p <- pp[ee]
    epsi <- Eps[ee]
    c1 <- log(1 + 1/p)
    c2 <- 2 * (p + epsi) * c1
    w0 <- (1/(epsi * c2)) * rep(1, m*q)
    flg <- 1

    while (1) {
      k <- k + 1

      # Acceleration double step
      V1 <- spEigenMMupdate(V, d, rho, Vx, sv2, w0, c1, p, epsi, m, q)
      V2 <- spEigenMMupdate(V1, d, rho, Vx, sv2, w0, c1, p, epsi, m, q)
      R <- V1 - V
      U <- V2 - V1 - R
      a <- min(-norm(abs(R), type = "F") / norm(abs(U), type = "F"), -1)

      # backtracking loop to ensure feasibility
      while (1) {
        V0 <- V - 2*a*R + a^2 * U

        # Projection
        s <- svd(V0)
        V0 <- s$u %*% h(s$v)

        g[abs(V0) <= epsi] <- abs(V0[abs(V0) <= epsi])^2 / (epsi*c2)
        g[abs(V0) > epsi] <- log((p + abs(V0[abs(V0) > epsi]))/(p + epsi))/c1 + epsi/c2
        if (data)
          F_v[k] <- colSums(abs(X %*% V0)^2) %*% d - colSums(g) %*% rho
        else
          F_v[k] <- Re(diag(h(V0) %*% X %*% V0)) %*% d - colSums(g) %*% rho

        if (flg == 0 && F_v[k] * (1 + sign(F_v[k])*1e-9) <= F_v[max(k - 1, 1)])
          a <- (a-1)/2
        else {
          V <- V0
          break
        }
      }

      # Stopping criterion
      if (flg == 0) {
        rel_change <- (abs(F_v[k] - F_v[k - 1]) / max(1, abs(F_v[k - 1]) ) ) # relative change in objective
        if (rel_change <= tol[ee] || k >= max_iter)
          break
      }
      flg <- 0
    }
  }

  V[abs(V) < thres] <- 0  # threshold
  nrm <- 1 / sqrt(colSums(abs(V)^2))
  V <- matrix(rep(nrm, m), ncol = q) * V

  return(list(vectors = V, standard_vectors = Vx[, 1:q], values = ifelse(rep(data,q), sv2[1:q] / (nrow(X) - 1), sv2[1:q])))
}


# MM update in ecah iteration
spEigenMMupdate <- function(V, d, rho, Vx, sv2, w0, c1, p, epsi, m, q) {

  # weights
  w <- w0
  absV <- abs(V)
  ind <- c(absV) > epsi
  w[ind] <- (0.5/c1) / (absV[ind]^2 + p*absV[ind])

  # calculation of H, G
  H <- (matrix(w, ncol=q) - rep(apply(matrix(w, ncol=q), 2, max), each = m)) * V * rep(rho, each = m)
  G <- Vx %*% ( (h(Vx) %*% (V * rep(d, each = m))) * matrix(rep(sv2, q), ncol = q) )  # Vx %*% diag(sv2) %*% t(Vx) %*% V * diag(sv2)

  # update
  s <- svd(G - H)

  return (Uk = s$u %*% h(s$v))
}


# Hermitian
h <- function(x) {
  return(Conj(t(x)))
}
