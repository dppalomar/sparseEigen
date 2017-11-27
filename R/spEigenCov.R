#' Covariance Matrix Estimation with Sparse Eigenvectors
#'
#' Computes joitly the sparse (orthogonal) eigenvectors and the eigenvalues of numeric covariance matrices.
#'
#' @param S m-by-m covariance matrix. It is required that \code{S} is full-rank. Both real and complex matrices are accepted.
#' @param q number of sparse eigenvectors.
#' @param rho sparsity weight factor. Any nonnegative number (suggested range [0,1]).
#' @param thres threshold value. All the entries of the sparse eigenvectors less or equal to \code{thres} are set to 0. The default value is \code{1e-9}.
#' @return A list with the following components:
#' \item{\code{vectors}  }{m-by-m matrix, columns corresponding to eigenvectors.}
#' \item{\code{values}  }{m-by-1 vector corresponding to eigenvalues.}
#' @author Konstantinos Benidis and Daniel P. Palomar
#' @references
#' K. Benidis, Y. Sun, P. Babu, D.P. Palomar "Orthogonal Sparse PCA and Covariance Estimation via Procrustes Reformulation,"
#' IEEE Transactions on Signal Processing, vol 64, no. 23, pp. 6211-6226, Dec. 2016.
#' @examples
#' library(sparseEigen)
#' n <- 600  # samples
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
#' # generate dats
#' X <- MASS::mvrnorm(n, rep(0, m), R)  # random data with underlying sparse structure
#'
#' # standard and spase estimation
#' res_standard <- eigen(cov(X))
#' res_sparse <- spEigenCov(cov(X), q, rho)
#'
#' # show inner product between estimated eigenvectors and originals (the closer to 1 the better)
#' abs(diag(t(res_standard$vectors) %*% V[, 1:q]))  #for standard estimated eigenvectors
#' abs(diag(t(res_sparse$vectors) %*% V[, 1:q]))    #for sparse estimated eigenvectors
#'
#' # show error between estimated and true covariance
#' norm(cov(X) - R, type = 'F') #for sample covariance matrix
#' norm(res_sparse$cov - R, type = 'F') #for covariance with sparse eigenvectors
#' @export
spEigenCov <- function(S, q = 1, rho = 0.5, thres = 1e-9) {
  max_iter <- 3000 # maximum MM iterations

  ######## error control  #########
  S <- as.matrix(S)
  m <- ncol(S)
  if (m == 1) stop("Data is univariate!")
  if (q > m) stop("The number of sparse eigenvectors q should not be larger than m.")
  if (anyNA(S) || anyNA(q) || anyNA(rho)) stop("This function cannot handle NAs.")
  if (!isSymmetric.matrix(S)) stop("The covariance matrix is not symmetric")
  if ( (q %% 1) != 0 || q <= 0) stop("The input argument q should be a positive integer.")
  if (rho <= 0) stop("The input argument rho should be positive.")
  #################################

  # EVD
  S_evd <- eigen(S)
  if (any(S_evd$values <= 0)) stop("The covariance matrix is not PSD.")
  if (sum(S_evd$values > 1e-9) < m) stop("The covariance matrix is low-rank.")
  V <- S_evd$vectors
  Xi <- S_evd$values
  Xi_max <- Xi[1]
  S_hat <- S - (Xi_max + 1e-6) * diag(rep(1, m))

  # Sparsity parameter rho
  rho <- rho * 10 * max(sqrt(sum(abs(S)^2))) %*% Xi[1:q]/Xi[1]

  # Preallocation
  F_v <- matrix(0, max_iter, 1) # objective value
  g <- matrix(0, m, q)

  # Decreasing epsilon, p
  K <- 5
  p1 <- 1 # first value of p
  pk <- 5 # last value of p
  gamma <- (pk / p1) ^ (1 / K)
  pp <- p1 * gamma ^ (0:K)
  pp <- 10 ^ (-pp)

  tol <- pp * 1e-1 # tolerance for convergence
  Eps <- pp * 1e-1 # epsilon

  ######################### MM LOOP #########################
  k <- 0  # MM iteration counter
  for (ee in 1:(K + 1)) {
    p <- pp[ee]
    epsi <- Eps[ee]
    c1 <- log(1 + 1/p)
    c2 <- 2 * (p + epsi) * c1
    w0 <- (1/(epsi * c2)) * rep(1, m*q)
    flg <- 1

    while (1) {
      k <- k + 1

      # Compute quantity A since it is used often
      A <- V %*% diag(1/Xi)

      # eigenvector update
      V <- eigvecUpdate(V, A, S_hat, rho, w0, c1, p, epsi, m, q)

      # eigenvalue update
      G <- -h(A) %*% S_hat %*% A  # this can be faster since we need only diag()
      Xi <- eigvalAlgo(Re(diag(G)), Xi_max, q)

      # Objective
      Vtmp <- abs(V[,1:q])
      g[abs(Vtmp) <= epsi] <- Vtmp[Vtmp <= epsi]^2 / (epsi*c2)
      g[abs(Vtmp) > epsi] <- log((p + Vtmp[Vtmp > epsi])/(p + epsi))/c1 + epsi/c2

      F_v[k] <- sum(log(Xi)) +  sum(diag(h(V) %*% S %*% V) * (1/Xi)) + rho %*% colSums(g)

      # Stopping criterion
      if (flg == 0) {
        rel_change <- abs(F_v[k] - F_v[k - 1]) / max(1, abs(F_v[k - 1])) # relative change in objective
        if ( (rel_change <= tol[ee]) | (k >= max_iter) ) {
          break
        }
      }
      flg <- 0
    }
  }

  V[abs(V) < thres] <- 0 # threshold
  nrm <- 1 / sqrt(colSums(abs(V)^2))
  V <- matrix(rep(nrm, m), ncol = m) * V

  return(list(vectors = V, values = Xi, cov = V %*% diag(Xi) %*% h(V)))
}


# Hermitian
h <- function(x) {
  return(Conj(t(x)))
}


# Eigenvector update
eigvecUpdate <- function(V, A, S_hat, rho, w0, c1, p, epsi, m, q) {
  w <- w0
  absV <- abs(V[,1:q])
  ind <- c(absV) > epsi
  w[ind] <- (0.5/c1) / (absV[ind]^2 + p*absV[ind])
  w <- c(w, rep(0, m*(m-q)))

  H1 <- (matrix(w, ncol = m) - rep(apply(matrix(w, ncol = m), 2, max), each = m)) * V * rep(c(rho,rep(0, m-q)), each = m)
  H2 <- S_hat %*% A

  # Procrustes update
  s <- svd(-(H1 + H2))

  return (Uk = s$u %*% h(s$v))
}


# Eigenvalue update algorithm
eigvalAlgo <- function(g, x, q) {
  m <- length(g)
  g_new <- g
  mult <- 0

  ##### Case 1: Many sparse eigenvectors (q>1) #####
  if (q > 1) {
    while (1) {
      flg <- 0
      i <- 1

      # Swaps among the first q-1.
      # If q is included in block swap then check the q+1:m last.
      while (i <= (q-1)) {
        if (g[i] >= g[i+1]) {

          if (i == (q-1)) {
            flg <- 1
          }

          # check block swaps
          j <- i + 1
          while (j <= (q-1)) {
            if (g[j] >= g[j+1]) {
              if (j == q-1) {
                flg <- 1
              }
              j <- j + 1
            }
            else {
              break
            }
          }

          # Swaps in the first q-1
          if (flg == 0) {
            g_new[i:j] <- 1/(j-i+1) * sum(g[i:j])
          }
          # Swaps in the first q-1 including the q-th. Check all the q+1:m
          else {
            swapInd <- which(g[q] >= g[(q+1):m]) + q

            if (length(swapInd) == 0) {
              g_new[i:j] <- 1/(j-i+1) * sum(g[i:j])
              break
            }

            redInd <- swapInd
            # Keep in the set the ones that are equal to a(q) from previous iterations
            if (mult == 1) {
              ind <- which(g[q] == g[(q+1):m]) + q
              ind <- c(seq(i:q), ind)
            }
            else {
              ind <- seq(i:q)
            }

            while (1) {
              if (length(redInd) > 0) { # to avoid warning
                indMin <- which(g[redInd] == min(g[redInd]))
                ind <- unique(c(ind, redInd[indMin]))
                redInd <- redInd[-indMin]
                tmp <- 1/length(ind) * sum(g[ind])
              }

              if (length(redInd) > 0) {
                if (sum(tmp > g[redInd]) == 0) {
                  g_new[ind] = 1/length(ind) * sum(g[ind])
                  break
                }
              }
              else {
                g_new[ind] = 1/length(ind) * sum(g[ind])
                break
              }
            }

            i <- j
          }
        }
        ##### Swaps only among the q-th and the (q+1:m) last #####
        else if ((i == (q-1)) & (flg == 0)) {
          swapInd <- which(g[q] >= g[(q+1):m]) + q

          redInd <- swapInd
          if (mult == 1) {
            ind <- which(g[q] == g[(q+1):m]) + q
            ind <- c(q, ind)
          }
          else {
            ind <- q
          }

          while (1) {
            if (length(redInd) > 0) {
              indMin <- which(g[redInd] == min(g[redInd]))
              ind <- unique(c(ind, redInd[indMin]))
              redInd <- redInd[-indMin]
              tmp <- 1/length(ind) * sum(g[ind])
            }

            if (length(redInd) > 0) {
              if (sum(tmp > g[redInd]) == 0) {
                g_new[ind] = 1/length(ind) * sum(g[ind])
                break
              }
            }
            else {
              g_new[ind] = 1/length(ind) * sum(g[ind])
              break
            }
          }
        }

        i <- i + 1
      }

      g <- g_new

      # Stopping criteria
      if (sum(order(g)[1:q] != seq(1:q)) == 0) {
        break
      }

      mult <- 1
    }

    phi <- (1 + sqrt(1 + 4 * x * g)) / (2 * x)
  }

  ##### Case 2: One sparse eigenvectors (q=1) #####
  if (q == 1) {
    while (1) {
      swapInd <- which(g[q] >= g[(q+1):m]) + q

      redInd <- swapInd
      if (mult == 1) {
        ind <- which(g[q] == g[(q+1):m]) + q
        ind <- c(q, ind)
      }
      else {
        ind <- q
      }

      while (1) {
        if (length(redInd) > 0) {
          indMin <- which(g[redInd] == min(g[redInd]))
          ind <- unique(c(ind, redInd[indMin]))
          redInd <- redInd[-indMin]
          tmp <- 1/length(ind) * sum(g[ind])
        }

        if (sum(tmp > g[redInd]) == 0) {
          g_new[ind] = 1/length(ind) * sum(g[ind])
          break
        }
        else {
          g_new[ind] = 1/length(ind) * sum(g[ind])
          break
        }
      }

      g <- g_new

      # Stopping criteria
      if (sum(order(g)[1:q] != seq(1:q)) == 0) {
        break
      }
    }

    phi <- (1 + sqrt(1 + 4 * x * g)) / (2 * x)
  }

  return (1 / phi)
}
