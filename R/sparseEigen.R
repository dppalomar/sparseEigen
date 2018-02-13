#' sparseEigen: Computation of Sparse Eigenvectors of a Matrix
#'
#'Computation of sparse eigenvectors of a matrix (aka sparse PCA)
#'with running time 2-3 orders of magnitude lower than existing methods and
#'better final performance in terms of recovery of sparsity pattern and
#'estimation of numerical values. Can handle covariance matrices as well as
#'data matrices with real or complex-valued entries. Different levels of
#'sparsity can be specified for each individual ordered eigenvector and the
#'method is robust in parameter selection. See vignette for a detailed
#'documentation and comparison, with several illustrative examples.
#'
#' @section Functions:
#' \code{\link{spEigen}}, \code{\link{spEigenCov}}
#'
#' @section Help:
#' For a quick help see the README \url{https://cran.r-project.org/web/packages/sparseEigen/README.html}.
#'
#' For more details see the vignette \url{https://cran.r-project.org/web/packages/sparseEigen/vignettes/SparseEigenvectors.pdf}.
#'
#' @author Konstantinos Benidis and Daniel P. Palomar
#'
#' @references
#' K. Benidis, Y. Sun, P. Babu, and D. P. Palomar, "Orthogonal Sparse PCA and Covariance Estimation via Procrustes Reformulation,"
#' \emph{IEEE Transactions on Signal Processing}, vol. 64, no. 23, pp. 6211-6226, Dec. 2016.
#'
#' @docType package
#' @name sparseEigen-package
NULL
