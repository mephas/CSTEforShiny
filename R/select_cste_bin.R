#' Select the optimal tuning parameters in CSTE estimation for binary outcome.  
#'
#' select lasso penalty parameter \eqn{\lambda} for \eqn{\beta_1} and
#'\eqn{\beta_2} in CSTE estimation.
#'  
#' 
#'@param x samples of covariates which is a \eqn{n*p} matrix.
#'@param y samples of binary outcome which is a \eqn{n*1} vector.
#'@param z samples of treatment indicator which is a \eqn{n*1} vector.
#'@param lam_seq a sequence for the choice of \eqn{\lambda}. 
#'@param beta_ini initial values for \eqn{(\beta_1', \beta_2')'}, default value is NULL.
#'@param nknots number of knots for the B-spline for estimating \eqn{g_1} and \eqn{g_2}.
#'@param max.iter maximum iteration for the algorithm.
#'@param eps numeric scalar \eqn{\geq} 0, the tolerance for the estimation 
#' of \eqn{\beta_1} and \eqn{\beta_2}. 
#'
#'@return A list which includes
#' \itemize{
#'    \item \code{optimal}: optimal cste within the given the sequence of \eqn{\lambda}.
#'    \item \code{bic}: BIC for the sequence of \eqn{\lambda}.
#'    \item \code{lam_seq}: the sequence of \eqn{\lambda} that is used.
#'    
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

select_cste_bin <- function(x, y, z, lam_seq, beta_ini = NULL, nknots = 1, max.iter = 2000, eps = 1e-3) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  out <- vector("list", length(lam_seq))
  beta1 <- matrix(0, p, length(lam_seq))
  beta2 <- matrix(0, p, length(lam_seq))
  for(i in 1:length(lam_seq)) {
    if(is.null(beta_ini)){
      beta_ini <- rep(normalize(rep(1, p)), 2)
    }
    if(i > 1) beta_ini <- c(out[[i-1]]$beta1, out[[i-1]]$beta2)
    out[[i]] <- cste_bin(x, y, z, lam = lam_seq[i], beta_ini = beta_ini, nknots = nknots, max.iter = max.iter, eps = eps)
    beta1[, i] <- out[[i]]$beta1
    beta2[, i] <- out[[i]]$beta2
    if(out[[i]]$flag | out[[i]]$df1 <= 2 | out[[i]]$df2 <= 2) {
      warnings("not all lambdas are used; decrease lambda")
      break
    }
  }
  df <- sapply(out[1:i], function(x) x$df)
  bic <- sapply(out[1:i], function(x) x$bic)
  loss <- sapply(out[1:i], function(x) x$loss)
  return(list(optimal = out[[which.min(bic)]], bic = bic, lam_seq = lam_seq[1:i], df = df, complete = out[1:i], beta1=beta1, beta2=beta2))
}