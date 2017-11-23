#' Sparse Spectral Decomposition of a Matrix
#'
#' Computes sparse (orthogonal) eigenvectors of covariance matrix or directly of data matrix.
#'
#' @param X m-by-m covariance matrix or n-by-m data matrix (n samples, m variables).
#' @param q number of eigenvectors to be estimated.
#' @param rho sparsity weight factor. Any nonnegative number (suggested range [0,1]).
#' @param data boolean variable. If \code{TRUE}, \code{X} is treated as a data matrix, else as a covariance matrix (default).
#' @param d vector with q weights. The default value is \code{rep(1, q)}.
#' @param V initial m-by-q matrix point. If not provided, the eigenvectors of the sample covariance matrix are used.
#' @param thres threshold value. All the entries of the sparse eigenvectors less or equal to \code{thres} are set to 0. The default value is \code{1e-9}.
#' @return A list with the following components:
#' \item{\code{vectors}}{m-by-q matrix, columns corresponding to the q leading sparse eigenvectors.}
#' \item{\code{standard_vectors}}{m-by-q matrix, columns corresponding to standard (non-sparse) leading eigenvectors.}
#' \item{\code{values}}{vector with the q leading eigenvalues in decreasing order.}
#' @author Konstantinos Benidis and Daniel P. Palomar
#' @references
#' K. Benidis, Y. Sun, P. Babu, D.P. Palomar "Orthogonal Sparse PCA and Covariance Estimation via Procrustes Reformulation,"
#' IEEE Transactions on Signal Processing, vol 64, no. 23, pp. 6211-6226, Dec. 2016.
#' @examples
#' library(sparseEigen)
#' n <- 100  # samples
#' m <- 500  # dimension
#' q <- 3  # number of sparse eigenvectors to be estimated
#' sp_card <- 0.1*m  # sparsity of each eigenvector
#'
#' # generate covariance matrix with sparse eigenvectors
#' V <- matrix(0, m, q)
#' V[cbind(seq(1, q*sp_card), rep(1:q, each = sp_card))] <- 1/sqrt(sp_card)
#' V <- cbind(V, matrix(rnorm(m*(m-q)), m, m-q))
#' V <- qr.Q(qr(V))  # orthogonalize eigenvectors
#' lmd <- c(100*seq(from = q, to = 1), rep(1, m-q))  # generate eigenvalues
#' R <- V %*% diag(lmd) %*% t(V)  # covariance matrix
#'
#' # generate data
#' X <- MASS::mvrnorm(n, rep(0, m), R)  # random data with underlying sparse structure
#'
#' # standardand sparse eigenvectors
#' res_standard <- eigen(cov(X))
#' res_sparse <- spEigen(cov(X), q, rho)
#'
#' # show inner product between estimated eigenvectors and originals (the closer to 1 the better)
#' abs(diag(t(res_standard$vectors) %*% V[, 1:q]))  #for standard estimated eigenvectors
#' abs(diag(t(res_sparse$vectors) %*% V[, 1:q]))    #for sparse estimated eigenvectors
#' @importFrom gmodels fast.svd
#' @export
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
  V_tld <- matrix(0, m, q)
  H <- matrix(0, m, q)
  F_v <- matrix(0, max_iter, 1)  # objective value
  g <- matrix(0, m, q)

  # Decreasing epsilon, p
  K <- 10
  p1 <- 1  # first value of p
  pk <- 7  # last value of p
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

      #-------------------------------------#
      # First iteration of the acceleration #

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
      s1 <- fast.svd(G - H)
      V1 <- s1$u %*% t(s1$v)

      #--------------------------------------#
      # Second iteration of the acceleration #

      # weights
      w <- w0
      ind <- abs(c(V1)) > epsi
      w[ind] <- (0.5 / c1) / (V1[ind] ^ 2 + p * abs(V1[ind]))

      # MM
      for (i in 1:q) {
        w_tmp <- w[( (i - 1) * m + 1):(i * m)]
        V_tld[, i] <- V1[, i] * d[i]
        H[, i] <- (w_tmp - max(w_tmp) * rep(1, m)) * V1[, i] * rho[i]
      }

      G <- svd_x$v %*% ( (t(svd_x$v) %*% V_tld) * matrix(rep(sv2, q), ncol = q) )

      # update
      s2 <- fast.svd(G - H)
      V2 <- s2$u %*% t(s2$v)

      #--------------#
      # Acceleration #

      R <- V1 - V
      U <- V2 - V1 - R
      a <- min(-norm(R, type = "F") / norm(U, type = "F"), -1)

      # backtracking loop
      while (1) {
        V0 <- V - 2 * a * R + a ^ 2 * U

        # Projection
        s3 <- fast.svd(V0)
        V0 <- s3$u %*% t(s3$v)

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
