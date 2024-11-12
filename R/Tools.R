#' tool functions
#' 
#' @param x samples of biomarker (or covariate) which is a \eqn{n*1} vector 
#' and should be scaled between 0 and 1.
#' @param x0 biomarkers at a fixed point.
#' @param h kernel bandwidth.
#' @param z treatment indicators which is \eqn{n*K} matrix.
#' @param input parameters with dimension (2*K+1).   
#' @param R indicator matrix of individuals at risk prior to \eqn{y_i}.
#' @param s censoring indicator which is a \eqn{n*1} vector.
#' @param sep parameter for trapezoid rule.
#' @param sloped the first derivative of \eqn{g}.
#' @param l contraction vector with dimension \eqn{K}.
#' @param beta0,beta0dot the value of beta and its first derivative at \eqn{x0}.
#' @param g0dot the first derivative of \eqn{g0} at x0. 
#' @param stdel,stgam,std store the values of \code{beta0}, \code{beta0dot},
#'  and \code{godot} at each \eqn{x_i}. 
#' @param m number of turns of resampling.  
#' @param i the \eqn{i}-th observation. 
#' @param alpha the (1-alpha)-confidence level of SCB. 
#' @param thres the output of \code{getthres} function. 
#' 
#' @return Intermediate results.
#' @name tool
#' @keywords internal
NULL

# transfer z to z*
#' @rdname tool
ztrans <- function(z, x, x0,h){
  zstar = cbind(z,diag((x-x0)/h) %*% z,(x-x0)/h)
  return(zstar)
}

#' @rdname tool
lpl <- function(input,x,x0,R,z,s,h){
  n = nrow(z)
  p = ncol(z)
  del = input[1:p]
  gam = input[(p+1):(2*p)]
  d = input[2*p +1]
  myker = gaussK((x-x0)/h)/h
  temp = z%*%del+diag(x-x0)%*%z%*%gam+d*(x-x0)
  return(-1/n*s%*%diag(myker)%*%(temp-log(R%*%diag(myker)%*%exp(temp))))
}


#' @rdname tool
g <- function(x0,sep,sloped){
  
  #trapezoid rule
  if(x0 == 1){
    return((sum(sloped)-0.5*sloped[1]-0.5*sloped[sep])/sep)
  }
  b = floor(x0*sep)+1
  temp = sum(sloped[1:b])-0.5*sloped[1]-0.5*sloped[b]
  pp = (x0-(b-1)*1/sep)*sep
  height = (1-pp)*sloped[b]+pp*sloped[b+1]
  temp = temp+(sloped[b]+height)*pp/2
  return(temp/sep)
}

#' @rdname tool
cntcov <- function(l,beta0,beta0dot,g0dot,x,x0,R,Z,s,h,sep,sloped){
  n = nrow(Z)
  p = ncol(Z)
  Zstar = ztrans(Z,x,x0,h)
  g0 = g(x0,sep,sloped)
  myker = gaussK((x-x0)/h)/h
  theta0 = c(beta0,h*beta0dot,h*g0dot)
  ns0 <- as.vector(R%*%diag(myker)%*%exp(Zstar%*%theta0+g0*rep(1,n)))
  ns1 <- R%*%diag(as.vector(myker*exp(Zstar%*%theta0+g0*rep(1,n))))%*%Zstar
  left <- 1/n*t(Zstar)%*%diag(as.vector((s*myker/ns0)%*%R*myker*exp(as.vector(Zstar%*%theta0+g0*rep(1,n)))))%*%Zstar
  right <- 1/n*t(ns1)%*%diag(s*myker/ns0/ns0)%*%ns1
  A <- left-right
  #B <- 1/n*(t(Zstar)%*%myker-t(ns1)%*%(myker/ns0))
  Pi <- 1/(n*h*sqrt(pi))*t(Zstar-diag(1/ns0)%*%ns1)%*%diag(s*myker*myker)%*%(Zstar-diag(1/ns0)%*%ns1) #checked
  mycov = 1/n/h*solve(A)%*%Pi%*%solve(A)
  return(as.numeric(l%*%mycov[1:p,1:p]%*%l))
}

#' @rdname tool
sampleQ <- function(i,l,stdel,stgam,std,x,R,Z,s,h,sep,sloped,m){
  n = nrow(Z)
  p = ncol(Z)
  lstar = c(l,rep(0,p+1))
  g0 = g(x[i],sep,sloped)
  Zstar = ztrans(Z,x,x[i],h)
  myker = gaussK((x-x[i])/h)/h
  G = matrix(rnorm(n*m),n,m)
  beta0 = as.vector(stdel[i,])
  beta0dot = as.vector(stgam[i,])
  g0dot = as.numeric(std[i])
  ns0 <- as.vector(R%*%diag(myker)%*%exp(Z%*%beta0+g0*rep(1,n)))
  ns1 <- R%*%diag(as.vector(myker*exp(Z%*%beta0+g0*rep(1,n))))%*%Zstar
  left <- 1/n*t(Zstar)%*%diag(as.vector((s*myker/ns0)%*%R*myker*exp(as.vector(Z%*%beta0+g0*rep(1,n)))))%*%Zstar
  right <- 1/n*t(ns1)%*%diag(s*myker/ns0/ns0)%*%ns1
  I <- left-right
  U <- 1/n*t(Zstar-diag(1/ns0)%*%ns1)%*%diag(myker*s)%*%G
  Q <- as.vector(lstar%*%solve(I)%*%U)
  mycov <- cntcov(l,beta0,beta0dot,g0dot,x,x[i],R,Z,s,h,sep,sloped)
  # return(abs(Q/sqrt(mycov)))
  return(abs(Q/sqrt(mycov)*sqrt(n*h)))
  
}

#' @rdname tool
getthres <- function(alpha,l,stdel,stgam,std,x,R,Z,s,h,sep,sloped,m){
  n = nrow(Z)
  tempfun <- function(t){return(sampleQ(t,l,stdel,stgam,std,x,R,Z,s,h,sep,sloped,m))}
  result <- matrix(sapply(1:n,tempfun),n,m,byrow=TRUE)
  Qbar <- apply(result,2,max)
  thres <- quantile(Qbar,1-alpha)
  return(thres)
}

#' @rdname tool
get.bound <- function(thres,l,stdel,stgam,std,x,R,Z,s,h,sep,sloped){
  n <- nrow(Z)
  p <- ncol(Z)
  tempfun <- function(i){
    beta0 = as.vector(stdel[i,])
    beta0dot = as.vector(stgam[i,])
    g0dot = as.numeric(std[i])
    return(cntcov(l,beta0,beta0dot,g0dot,x,x[i],R,Z,s,h,sep,sloped))}
  covs <- sapply(1:n,tempfun)
  return(cbind(stdel%*%l-thres*sqrt(covs)/sqrt(n*h),stdel%*%l,stdel%*%l+thres*sqrt(covs)/sqrt(n*h)))
}
