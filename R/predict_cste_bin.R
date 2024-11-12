#' Predict the CSTE curve of new data for binary outcome. 
#'
#' Predict the CSTE curve of new data for binary outcome.
#' 
#' 
#'@param obj a S3 class of cste.
#'@param newx samples of covariates which is a \eqn{m*p} matrix.
#'
#'@return A S3 class of cste which includes 
#' \itemize{
#'    \item \code{g1}: predicted \eqn{g_1(X\beta_1)}. 
#'    \item \code{g2}: predicted \eqn{g_2(X\beta_2)}.
#'    \item \code{B1}: the B-spline basis for estimating \eqn{g_1}.
#'    \item \code{B2}: the B-spline basis for estimating \eqn{g_2}. 
#' }
#'
#' @references
#' Guo W., Zhou X. and Ma S. (2021).
#' Estimation of Optimal Individualized Treatment Rules
#' Using a Covariate-Specific Treatment Effect Curve with 
#' High-dimensional Covariates,
#' \emph{Journal of the American Statistical Association}, 116(533), 309-321
#' 
#' @seealso \code{\link{cste_bin}}

predict_cste_bin <- function(obj, newx) {
  # type <- match.arg(type)
  if(missing(newx)) {
    out <- obj$B1 %*% obj$delta1
    newB <- NULL
  } else {
    u1 <- pu(newx, obj$beta1)
    u2 <- pu(newx, obj$beta2)
    eta1 <- u1$u
    newB <- bsplineS(eta1, breaks = quantile(eta1, obj$knots))
    g1 <- newB %*% obj$delta1
    eta2 <- u2$u
    newB2 <- bsplineS(eta2, breaks = quantile(eta2, obj$knots))
    g2 <- newB2 %*% obj$delta2
  }
  return(list(g1 = g1, g2 = g2, B1 = newB, B2 = newB2))
}