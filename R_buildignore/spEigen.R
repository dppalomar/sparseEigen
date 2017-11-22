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
    d <- seq(from = 1, to = 0.1, length.out = q)

  # Sparsity parameter rho
  svd_x <- fast.svd(X)
  if (data) {
    sv2 <- svd_x$d ^ 2
    rho <- rho * max(colSums(X ^ 2)) * (sv2[1:q] / sv2[1]) * d
  }
  else {
    sv2 <- svd_x$d
    rho <- rho * max(diag(X)) * (sv2[1:q] / sv2[1]) * d
  }

  # Input parameter V: initial point
  if (is.na(V))
    V <- svd_x$v[, 1:q]

  # Preallocation
  F_v <- matrix(0, max_iter, 1)  # objective value
  g <- matrix(0, m, q)

  # Decreasing epsilon, p
  K <- 10
  p1 <- 1  # first value of -log(p)
  pk <- 7  # last value of -log(p)
  gamma <- (pk / p1) ^ (1 / K)
  pp <- p1 * gamma ^ (0:K)
  pp <- 10 ^ (-pp)

  tol <- pp * 1e-2  # tolerance for convergence
  Eps <- pp  # epsilon


  ######################### MM LOOP #########################
  k <- 0  # MM iteration counter
  for (ee in 1:(K + 1)) {  # loop for approximation based on p & epsilon
    p <- pp[ee]
    epsi <- Eps[ee]
    c1 <- log(1 + 1 / p)
    c2 <- 2 * (p + epsi) * c1
    w0 <- (1 / (epsi * c2)) * rep(1, m * q)
    flg <- 1

    while (1) {
      k <- k + 1

      # First iteration of the acceleration
      V1 <- spEigenMMupdate(V, d, rho, svd_x, sv2, w0, c1, p, epsi, m, q)

      # Second iteration of the acceleration
      V2 <- spEigenMMupdate(V1, d, rho, svd_x, sv2, w0, c1, p, epsi, m, q)

      #--------------#
      # Acceleration #

      R <- V1 - V
      U <- V2 - V1 - R
      a <- min(-norm(R, type = "F") / norm(U, type = "F"), -1)

      # backtracking loop
      while (1) {
        V0 <- V - 2 * a * R + a ^ 2 * U

        # Projection
        s <- fast.svd(V0)
        V0 <- s$u %*% t(s$v)

        g[abs(V0) <= epsi] <- V0[abs(V0) <= epsi] ^ 2 / (epsi * c2)
        g[abs(V0) > epsi] <- (log( (p + abs(V0[abs(V0) > epsi]) ) / (p + epsi) )
                              / c1 + epsi / c2)
        if (data)
          F_v[k] <- colSums( (X %*% V0) ^ 2) %*% d - colSums(g) %*% rho
        else
          F_v[k] <- diag(t(V0) %*% X %*% V0) %*% d - colSums(g) %*% rho

        if (flg == 0 && F_v[k] * (1 + sign(F_v[k]) * 1e-9) <= F_v[max(k - 1, 1)])
          a <- (a - 1) / 2
        else {
          V <- V0
          break
        }
      }

      # Stopping criterion
      if (flg == 0) {
        rel_change <- (abs(F_v[k] - F_v[k - 1])
                       / max(1, abs(F_v[k - 1]) ) ) # relative change in objective
        if (rel_change <= tol[ee] || k >= max_iter)
          break
      }
      flg <- 0
    }
  }

  V[abs(V) < thres] <- 0  # threshold
  nrm <- 1 / sqrt(colSums(V^2))
  V <- matrix(rep(nrm, m), ncol = q) * V

  return(list(vectors = V, standard_vectors = svd_x$v[, 1:q], values = ifelse(rep(data,q), sv2[1:q] / (nrow(X) - 1), sv2[1:q])))
}



spEigenMMupdate <- function(V, d, rho, svd_x, sv2, w0, c1, p, epsi, m, q) {
  V_tld <- matrix(0, m, q)
  H <- matrix(0, m, q)

  # weights
  w <- w0
  ind <- abs(c(V)) > epsi
  w[ind] <- (0.5 / c1) / (V[ind] ^ 2 + p * abs(V[ind]))

  # MM
  for (i in 1:q) {
    w_tmp <- w[( (i - 1) * m + 1):(i * m)]
    V_tld[, i] <- V[, i] * d[i]
    H[, i] <- (w_tmp - max(w_tmp) * rep(1, m)) * V[, i] * rho[i]
  }

  G <- svd_x$v %*% ( (t(svd_x$v) %*% V_tld) * matrix(rep(sv2, q), ncol = q) )

  # update
  s <- fast.svd(G - H)

  return (s$u %*% t(s$v))
}
