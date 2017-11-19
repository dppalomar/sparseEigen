#' Covariance Matrix Estimation with Sparse Eigenvectors
#'
#' Computes joitly the sparse (orthogonal) eigenvectors and the eigenvalues of numeric covariance matrices.
#'
#' @param S m-by-m covariance matrix.
#' @param q number of sparse eigenvectors.
#' @param rho sparsity weight factor. Any nonnegative number. Suggested range [0, 1].
#' @param thres threshold value. All the entries of the sparse eigenvectors less or equal to \code{thres} are set to 0. The default value is \code{1e-9}.
#' @return A list with the following components:
#' \item{\code{vectors}  }{m-by-m matrix, columns corresponding to eigenvectors.}
#' \item{\code{values}  }{m-by-1 vector corresponding to eigenvalues.}
#' @author Konstantinos Benidis and Daniel P. Palomar
#' @references
#' K. Benidis, Y. Sun, P. Babu, D.P. Palomar "Orthogonal Sparse PCA and Covariance Estimation via Procrustes Reformulation,"
#' IEEE Transactions on Signal Processing, vol 64, no. 23, pp. 6211-6226, Dec. 2016.
#' @examples
#' @export
spEigenCov <- function(S, q = 1, rho = 0.5, thres = 1e-9) {

  m <- ncol(S)

  ######## error control  #########
  if (m == 1) stop("Data is univariate!")
  if (qr(S)$rank < m) stop("The covariance is low-rank.")
  if (q > m) stop("The number of sparse eigenvectors q should not be larger than m.")
  if (anyNA(S) || anyNA(q) || anyNA(rho)) stop("This function cannot handle NAs.")
  if ( (q %% 1) != 0 || q <= 0) stop("The input argument q should be a positive integer.")
  if (rho <= 0) stop("The input argument rho should be positive.")
  #################################

  # EVD
  S_evd <- eigen(S)
  V <- S_evd$vectors
  Xi <- S_evd$values
  Xi_max <- Xi[1]
  S_hat <- S - (Xi_max + 1e-6) * diag(rep(1, m))

  # rho
  rho <- rho * 10 * max(sqrt(sum(S ^ 2))) %*% Xi[1:q] / Xi[1]

  # MM Parameters
  k <- 0 # MM iteration counter
  max_iter <- 3000 # maximum MM iterations

  # Preallocation
  F_v <- matrix(0, max_iter, 1) # objective value
  H1 <- matrix(0, m, m)
  H2 <- matrix(0, m, m)
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


  for (ee in 1:(K + 1)) {
    p <- pp[ee]
    epsi <- Eps[ee]
    c1 <- log(1 + 1 / p)
    c2 <- 2 * (p + epsi) * c1
    w0 <- (1 / (epsi * c2)) * rep(1, m * q)
    flg <- 1

    while (1) {
      k <- k + 1

      # Compute quantity A since it is used often
      A <- V %*% diag(1/Xi)

      ##### EIGENVECTOR UPDATE #####
      # weights
      w <- w0
      ind <- abs(c(V[,1:q])) > epsi
      Vtmp <- V[,1:q]
      w[ind] <- (0.5 / c1) / (Vtmp[ind] ^ 2 + p * abs(Vtmp[ind]))
      w <- c(w, rep(0, m*(m-q)))

      for (i in 1:q) {
        w_tmp <- w[( (i - 1) * m + 1):(i * m)]
        H1[, i] <- (w_tmp - max(w_tmp) * rep(1, m)) * V[, i] * rho[i]
      }

      H2 <- S_hat %*% A

      # Procrustes update
      s <- fast.svd(-(H1 + H2))
      V <- s$u %*% t(s$v)

      ##### EIGENVALUE UPDATE #####
      G <- -t(A) %*% S_hat %*% A  # this can be faster since we need only diag()
      Xi <- eigvalAlgo(diag(G), Xi_max, q)

      # Objective
      Vtmp <- V[,1:q]
      g[abs(Vtmp) <= epsi] <- Vtmp[abs(Vtmp) <= epsi] ^ 2 / (epsi * c2)
      g[abs(Vtmp) > epsi] <- (log( (p + abs(Vtmp[abs(Vtmp) > epsi]) ) / (p + epsi) )
                              / c1 + epsi / c2)

      F_v[k] <- (sum(log(Xi)) +  sum(diag(t(V) %*% S %*% V) * (1/Xi))
                 + rho %*% colSums(g))

      # Stopping criterion
      if (flg == 0) {
        rel_change <- (abs(F_v[k] - F_v[k - 1])
                       / max(1, abs(F_v[k - 1]) ) ) # relative change in objective
        if ( (rel_change <= tol[ee]) | (k >= max_iter) ) {
          break
        }
      }
      flg <- 0
    }
  }

  V[abs(V) < thres] <- 0 # threshold
  nrm <- 1 / sqrt(colSums(V ^ 2))
  V <- matrix(rep(nrm, m), ncol = m) * V

  return(list(vectors = V, values = Xi))
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
