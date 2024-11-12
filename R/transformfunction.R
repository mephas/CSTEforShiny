#' tool functions
#' 
#' @param x numeric vector or matrix.
#' @param y,u,u1,u2,beta,atrisk numeric vector.
#' @param fit a S3 class of cste.
#' @param off numeric value, offset.
#' @param lam numeric value, penalty parameter.
#' @param pen hyper-parameter that used in MCP and SCAD penalty functions.
#' 
#' @return Intermediate results.
#' @name tool
#' @keywords internal
NULL

# u transformation
#' @rdname tool
pu <- function(x, beta) {
  x <- as.matrix(x)          
  d <- (length(beta) + 1)/2
  v <- x%*%beta
  # alpha <- quantile(apply(x, 1, function(x) sqrt(sum(x^2))), 0.95)
  alpha <- quantile( sqrt( rowSums(x*x) ), 0.95) 
  t <- (v + alpha)/(2 * alpha)
  u <- pbeta(t, d, d)
  deriv <- dbeta(t, d, d)/(2*alpha)
  return(list(u = u, deriv = drop(deriv)))
}

#' @rdname tool
pu_inv <- function(x, beta, u) {
  d <- (length(beta) + 1)/2
  # alpha <- quantile(apply(x, 1, function(x) sqrt(sum(x^2))), 0.95)
  alpha <- quantile( sqrt( rowSums(x*x)), 0.95) 
  t <- qbeta(u, d, d)
  return(2 * alpha * t - alpha)
}

# normalize vector
#' @rdname tool
normalize <- function(x) x/sqrt(sum(x^2))

# logit generator
#' @rdname tool
logitinv <- function(x) exp(x)/(1+exp(x))

# sparse logistic regression without intercept
#' @rdname tool
my_logit <- function(x, y, off = NULL, beta = NULL, lam = 0, pen = 2) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  if(is.null(off)) off <- rep(0, n)	
  if(is.null(beta)) beta <- rep(0, p)
  return(penC(x, y, off, beta, lam, pen))
}


# intermediate functions for survival analysis 
#' @rdname tool
my_surv <- function(x, y, atrisk, off = NULL, beta = NULL, lam = 0) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  if(is.null(off)) off <- rep(0, n)	
  if(is.null(beta)) beta <- normalize(rep(1, p))
  penal <- function(bt){
    loglik = sum(off) + sum(x%*%bt) - sum(log(atrisk%*%exp(off+x%*%bt)))
    return(-loglik + lam * sum(abs(bt)))
  }
  beta_new <- optim(par=beta, fn=penal)$par
  return(beta_new)
}

# intermediate functions for cste_bin
#' @rdname tool
prev_fit_cste <- function(u1, u2, fit) {
  # initialization
  z <- fit$z
  id0 <- (z == 0)
  id1 <- (z == 1)
  u1 <- u1[id1]
  u2 <- u2[id1]
  y <- (fit$y)[id1]
  n <- sum(id1)
  fit_g1 <- (fit$B1 %*% fit$delta1)[id1]
  fit_g2 <- (fit$B2 %*% fit$delta2)[id1]
  fit_g <- fit$g[id1]
  
  # number of interior knots
  nk <- max(1, floor(min(n^(0.25)*log(n)+1, n/8-0.5-1)))
  
  # calculate optimal bandwidth
  mx.dt = na.omit(cbind(u1, fit_g1))
  mx <- smooth.spline(mx.dt[,1], mx.dt[,2])
  mx_derv <- predict(mx, u1, deriv = 1)$y
  mx_derv_2 <- predict(mx, u1, deriv = 2)$y
  
  # estimate sigma_b^2(x)
  sigma_b <- exp(fit_g) / (1 + exp(fit_g))^2
  fit_sigma_b.dt = na.omit(cbind(u1, sigma_b))
  fit_sigma_b <- smooth.spline(fit_sigma_b.dt[,1], fit_sigma_b.dt[,2])
  
  # estimate sigma^2(x)
  sigma_x <- (y - logitinv(fit_g))^2
  fit_sigma_x.dt = na.omit(cbind(u1, sigma_x))
  fit_sigma_x <- smooth.spline(fit_sigma_x.dt[,1], fit_sigma_x.dt[,2])
  
  # estimate kernel density
  fit_f_x <- density(u1)
  f_x <- approx(fit_f_x$x, fit_f_x$y, u1)$y
  # h_x <- bw.nrd0(u1)
  # f_x <- sapply(u1, function(xx) mean(dnorm((xx - u1)/h_x)/h_x))
  
  # calculate bias
  Ck_d <- 1/(2*sqrt(pi))
  Ck_n <- 1/(4*sqrt(pi))
  
  D <- predict(fit_sigma_b, u1)$y * f_x
  v_sq <- Ck_d * f_x * predict(fit_sigma_x, u1)$y
  
  # v_sq[v_sq < 0] <- 0.001
  # D[D < 0] <- 0.001
  
  b3 <- exp(fit_g) * (exp(fit_g) - 1) / (exp(fit_g) + 1)^3
  fit_b3.dt = na.omit(cbind(u1, b3))
  fit_b3 <- smooth.spline(fit_b3.dt[,1], fit_b3.dt[,2])
  
  mu2k <- 1/(2*sqrt(2))
  bias_x <- mu2k * (mx_derv_2 * D + mx_derv * f_x * predict(fit_sigma_b, u1, deriv=1)$y - mx_derv^2 * f_x * predict(fit_b3, u1)$y)
  h_opt <- (mean((1/D) * v_sq * (1/D)) / (4 * sum((bias_x / D)^2)))^(0.2)
  
  return(list(h = h_opt, fit_sigma_x = fit_sigma_x, fit_sigma_b = fit_sigma_b, fit_f_x = fit_f_x, n = n, y = y, u1 = u1, u2 = u2, fit_g1 = fit_g1, fit_g2 = fit_g2))
}
