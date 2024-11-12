#' Estimate the CSTE curve for time to event outcome with right censoring.  
#'
#' Estimate the CSTE curve for time to event outcome with right censoring. 
#' The working model 
#' is \deqn{\lambda(t| X, Z) = \lambda_0(t) \exp(\beta^T(X)Z + g(X)),}
#' which implies that \eqn{CSTE(x) = \beta(x)}.   
#' 
#' 
#'@param x samples of biomarker (or covariate) which is a \eqn{n*1} vector 
#' and should be scaled between 0 and 1.
#'@param y samples of time to event which is a \eqn{n*1} vector.
#'@param z samples of treatment indicator which is a \eqn{n*K} matrix.
#'@param s samples of censoring indicator which is a \eqn{n*1} vector. 
#'@param h kernel bandwidth.
#'@return A \eqn{n*K} matrix, estimation of \eqn{\beta(x)}. 
#' @references
#' Ma Y. and Zhou X. (2017). 
#' Treatment selection in a randomized clinical trial via covariate-specific 
#' treatment effect curves, \emph{Statistical Methods in Medical Research}, 26(1), 124-141.
#' 
#' @seealso \code{\link{cste_surv_SCB}}

cste_surv <- function(x,y,z,s,h){
  n <- nrow(z)
  p <- ncol(z)
  sep <- 20
  myfun <- function(w){return(as.numeric(y>=w))}
  R <- matrix(sapply(y,myfun),n,n,byrow=TRUE)
  #preprocessing
  stdel <- matrix(0,nrow=n,ncol=p)
  for(i in 1:n){
    tempfun <- function(t){return(lpl(t,x,x[i],R,z,s,h))}
    # ans = optim(rep(0,2*p+1),tempfun)$par
    ans = nmk(rep(0,2*p+1),tempfun)$par
    stdel[i,] <- ans[1:p]
  }
  return(stdel)
}


