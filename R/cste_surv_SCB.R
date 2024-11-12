#' Calculate simultaneous confidence bands (SCB) of CSTE curve for time to event outcome with right censoring.
#' 
#' This function calculates simultaneous confidence bands of CSTE curve for time to event outcome with right censoring.
#' 
#' @param l contraction vector with dimension \eqn{K}.
#' @param x samples of biomarker (or covariate) which is a \eqn{n*1} vector 
#' and should be scaled between 0 and 1.
#' @param y samples of time to event which is a \eqn{n*1} vector.
#' @param z samples of treatment indicator which is a \eqn{n*K} matrix.
#' @param s samples of censoring indicator which is a \eqn{n*1} vector. 
#' @param h kernel bandwidth.
#' @param m number of turns of resampling.  
#' @param alpha the \eqn{(1-\alpha)}-confidence level of SCB. 
#'
#'@return A \eqn{n*3} matrix, estimation of \eqn{l^T \beta(x)} and its simultaneous confidence bands. 
#' @references
#' Ma Y. and Zhou X. (2017). 
#' Treatment selection in a randomized clinical trial via covariate-specific 
#' treatment effect curves, \emph{Statistical Methods in Medical Research}, 26(1), 124-141.
#' 
#' @seealso \code{\link{cste_surv}}

cste_surv_SCB <- function(l,x,y,z,s,h,m, alpha= 0.05){
  n <- nrow(z)
  p <- ncol(z)
  sep <- 20
  myfun <- function(w){return(as.numeric(y>=w))}
  R <- matrix(sapply(y,myfun),n,n,byrow=TRUE)
  #preprocessing
  stdel <- matrix(0,nrow=n,ncol=p)
  stgam <- matrix(0,nrow=n,ncol=p)
  std <- rep(0,n)
  for(i in 1:n){
    tempfun <- function(t){return(lpl(t,x,x[i],R,z,s,h))}
    # ans = optim(rep(0,2*p+1),tempfun)$par
    ans = nmk(rep(0,2*p+1),tempfun)$par
    stdel[i,] <- ans[1:p]
    stgam[i,] <- ans[(p+1):(2*p)]
    std[i] <- ans[2*p+1]
  }
  sloped <- rep(0,sep)
  for(i in 1:(sep+1)){
    tempfun <- function(t){return(lpl(t,x,1/sep*(i-1),R,z,s,h))}
    # ans = optim(rep(0,2*p+1),tempfun)$par
    ans = nmk(rep(0,2*p+1),tempfun)$par
    sloped[i] <- ans[2*p+1]
  }
  
  # derive quantile
  thres <- getthres(alpha,l,stdel,stgam,std,x,R,z,s,h,sep,sloped,m)
  return(get.bound(thres,l,stdel,stgam,std,x,R,z,s,h,sep,sloped))
}
