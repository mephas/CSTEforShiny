#' Estimate the CSTE curve for binary outcome.  
#'
#' Estimate covariate-specific treatment effect (CSTE) curve. Input data
#' contains covariates \eqn{X}, treatment assignment \eqn{Z} and binary outcome
#' \eqn{Y}. The working model is \deqn{logit(\mu(X, Z)) = g_1(X\beta_1)Z + g_2(X\beta_2),} 
#' where \eqn{\mu(X, Z) = E(Y|X, Z)}. The model implies that \eqn{CSTE(x) = g_1(x\beta_1)}.   
#' 
#'@param x samples of covariates which is a \eqn{n*p} matrix.
#'@param y samples of binary outcome which is a \eqn{n*1} vector.
#'@param z samples of treatment indicator which is a \eqn{n*1} vector.
#'@param beta_ini initial values for \eqn{(\beta_1', \beta_2')'}, default value is NULL.
#'@param lam value of the lasso penalty parameter \eqn{\lambda} for \eqn{\beta_1} and
#'\eqn{\beta_2}, default value is 0.
#'@param nknots number of knots for the B-spline for estimating \eqn{g_1} and \eqn{g_2}.
#'@param max.iter maximum iteration for the algorithm.
#'@param eps numeric scalar \eqn{\geq} 0, the tolerance for the estimation 
#' of \eqn{\beta_1} and \eqn{\beta_2}.  
#'
#'@return A S3 class of cste, which includes: 
#' \itemize{
#'    \item \code{beta1}: estimate of \eqn{\beta_1}.
#'    \item \code{beta2}: estimate of \eqn{\beta_2}.
#'    \item \code{B1}: the B-spline basis for estimating \eqn{g_1}.
#'    \item \code{B2}: the B-spline basis for estimating \eqn{g_2}.
#'    \item \code{delta1}: the coefficient of B-spline for estimating \eqn{g_1}.
#'    \item \code{delta2}: the coefficient for B-spline for estimating \eqn{g_2}.
#'    \item \code{iter}: number of iteration.
#'    \item \code{g1}: the estimate of \eqn{g_1(X\beta_1)}.
#'    \item \code{g2}: the estimate of \eqn{g_2(X\beta_2)}. 
#' }
#'@examples
#' ## Quick example for the cste
#'
#' library(mvtnorm)
#' library(sigmoid)
#' 
#' # --------  Example 1: p = 20 ---------  #  
#' ## generate data 
#' n <- 2000
#' p <- 20
#' set.seed(100)
#' 
#' # generate X
#' sigma <- outer(1:p, 1:p, function(i, j){ 2^(-abs(i-j)) } )
#' X <- rmvnorm(n, mean = rep(0,p), sigma = sigma)
#' X <- relu(X + 2) - 2
#' X <- 2 - relu(2 - X)
#' 
#' # generate Z
#' Z <- rbinom(n, 1, 0.5)
#' 
#' # generate Y
#' beta1 <- rep(0, p)
#' beta1[1:3] <- rep(1/sqrt(3), 3)
#' beta2 <- rep(0, p)
#' beta2[1:2] <- c(1, -2)/sqrt(5)
#' mu1 <- X %*% beta1
#' mu2 <- X %*% beta2
#' g1 <- mu1*(1 - mu1)
#' g2 <- exp(mu2)      
#' prob <- sigmoid(g1*Z + g2)
#' Y <- rbinom(n, 1, prob)
#' 
#' ## estimate the CSTE curve
#' fit <- cste_bin(X, Y, Z)
#' 
#' ## plot 
#' plot(mu1, g1, cex = 0.5, xlim = c(-2,2), ylim = c(-8, 3), 
#'      xlab = expression(X*beta), ylab = expression(g1(X*beta)))
#'      ord <- order(mu1)
#'      points(mu1[ord], fit$g1[ord], col = 'blue', cex = 0.5)
#'      
#' ## compute 95% simultaneous confidence band (SCB)
#' res <- cste_bin_SCB(X, fit, alpha = 0.05)
#' 
#' ## plot 
#' plot(res$or_x, res$fit_x, col = 'red', 
#'      type="l", lwd=2, lty = 3, ylim = c(-10,8),
#'      ylab=expression(g1(X*beta)), xlab = expression(X*beta), 
#'      main="Confidence Band")
#' lines(res$or_x, res$lower_bound, lwd=2.5, col = 'purple', lty=2)
#' lines(res$or_x, res$upper_bound, lwd=2.5, col = 'purple', lty=2)
#' abline(h=0, cex = 0.2, lty = 2)
#' legend("topleft", legend=c("Estimates", "SCB"), 
#'         lwd=c(2, 2.5), lty=c(3,2), col=c('red', 'purple'))
#'         
#'         
#' # --------  Example 2: p = 1 ---------  #  
#' 
#' ## generate data 
#' set.seed(15)
#' p <- 1
#' n <- 2000
#' X <- runif(n)
#' Z <- rbinom(n, 1, 0.5)
#' g1 <- 2 * sin(5*X) 
#' g2 <- exp(X-3) * 2
#' prob <- sigmoid( Z*g1 + g2)
#' Y <- rbinom(n, 1, prob)
#' 
#' ## estimate the CSTE curve
#' fit <- cste_bin(X, Y, Z)  
#' 
#' ## simultaneous confidence band (SCB)
#' X <- as.matrix(X)
#' res <- cste_bin_SCB(X, fit)  
#' 
#' ## plot 
#' plot(res$or_x, res$fit_x, col = 'red', type="l", lwd=2, 
#'      lty = 3, xlim = c(0, 1), ylim = c(-4, 4), 
#'      ylab=expression(g1(X)), xlab = expression(X), 
#'      main="Confidence Band")
#' lines(res$or_x, res$lower_bound, lwd=2.5, col = 'purple', lty=2)
#' lines(res$or_x, res$upper_bound, lwd=2.5, col = 'purple', lty=2)
#' abline(h=0, cex = 0.2)
#' lines(X[order(X)], g1[order(X)], col = 'blue', lwd = 1.5)
#' legend("topright", legend=c("Estimates", "SCB",'True CSTE Curve'), 
#' lwd=c(2, 2.5, 1.5), lty=c(3,2,1), col=c('red', 'purple','blue'))
#' 
#' @references
#' Guo W., Zhou X. and Ma S. (2021).
#' Estimation of Optimal Individualized Treatment Rules
#' Using a Covariate-Specific Treatment Effect Curve with 
#' High-dimensional Covariates,
#' \emph{Journal of the American Statistical Association}, 116(533), 309-321
#' 
#' @seealso \code{\link{cste_bin_SCB}, \link{predict_cste_bin}, \link{select_cste_bin}}


# cste estimation for binary outcome
cste_bin <- function(x, y, z, beta_ini = NULL, lam = 0, nknots = 1, max.iter = 200, eps = 1e-3) {
  x <- as.matrix(x)
  n <- dim(x)[1]
  p <- dim(x)[2]
  if(p==1) {
    B1 <- B2 <- bs(x, df = nknots+4, intercept = TRUE) 
    B <- cbind(z * B1, B2)
    # fit <- glm(y~0+B, family="binomial")
    fit <- glm.fit(B, y, family=binomial(link = 'logit'))  
    delta1 <- coef(fit)[1:(nknots+4)]
    delta2 <- coef(fit)[(nknots+5):(2*nknots+8)]
    geta <- z * B1 %*% delta1 + B2 %*% delta2
    g1 <- B1 %*% delta1
    g2 <- B2 %*% delta2
    loss <- -sum(log(1 + exp(geta))) + sum(y * geta)
    bic <- -2 * loss + (nknots+4) * log(n) 
    aic <- -2 * loss + 2 * (nknots+4) 
    out <- list(beta1 = 1, beta2 = 1, B1 = B1, B2 = B2, delta2 = delta2, delta1 = delta1, g = geta, x = x, y = y, z = z, nknots = nknots, p=p, g1=g1, g2=g2, bic=bic, aic=aic)
    class(out) <- "cste"
    return(out)
  } else {
    flag <- FALSE
    # if(is.null(truth)) 
    truth <- rep(p,2)
    if(is.null(beta_ini)) beta_ini <- c(normalize(rep(1, truth[1])), normalize(rep(1, truth[2])))
    else beta_ini <- c(beta_ini[1:truth[1]], beta_ini[(truth[1]+1):sum(truth)])
    beta_curr <- beta_ini
    conv <- FALSE
    iter <- 0
    # knots location
    knots <- seq(0, 1, length = nknots + 2)
    len.delta <- length(knots) + 2
    while(conv == FALSE & iter < max.iter) {
      iter <- iter + 1
      # step 1 fix beta to estimate g
      beta1 <- beta_curr[1:truth[1]]
      beta2 <- beta_curr[(truth[1]+1):length(beta_curr)]
      # u transformation
      u1 <- pu(x[,1:truth[1]], beta1)
      u2 <- pu(x[,1:truth[2]], beta2)
      eta1 <- u1$u
      eta2 <- u2$u
      # calculate B-spline basis
      B1 <- bsplineS(eta1, breaks = quantile(eta1, knots))
      B2 <- bsplineS(eta2, breaks = quantile(eta2, knots))
      B <- cbind(z*B1, B2)
      # estimate g1 and g2
      # fit.delta <- glm(y~0+B, family="binomial")
      fit.delta <- glm.fit(B, y, family=binomial(link = 'logit'))  
      delta <- drop(coef(fit.delta))
      delta[is.na(delta)] <- 0
      delta1 <- delta[1 : len.delta]
      delta2 <- delta[(len.delta + 1) : (2*len.delta)]
      # step 2 fix g to estimate beta
      # calculate first derivative of B-spline basis
      B_deriv_1 <- bsplineS(eta1, breaks = quantile(eta1, knots), nderiv = 1)
      B_deriv_2 <- bsplineS(eta2, breaks = quantile(eta2, knots), nderiv = 1)
      # calculate input 
      newx_1 <- z * drop((B_deriv_1*(u1$deriv))%*%delta1)*x[,1:truth[1]]
      newx_2 <- drop((B_deriv_2*(u2$deriv))%*%delta2)*x[,1:truth[2]]
      newx <- cbind(newx_1, newx_2)
      # calculate offset
      off_1 <- z * B1%*%delta1 - newx_1 %*% beta1
      off_2 <- B2%*%delta2 - newx_2 %*% beta2
      off <- off_1 + off_2
      # estimate beta
      beta <- my_logit(newx, y, off, lam = lam)
      if(sum(is.na(beta)) > 1){
        break
        stop("only 1 variable in betas; decrease lambda")
      }
      beta1 <- beta[1:truth[1]]
      beta2 <- beta[(truth[1]+1):length(beta_curr)]
      check <- c(sum(beta1!=0),sum(beta2!=0))
      if(min(check) <= 1) {
        stop("0 beta occurs; decrease lambda")
        flag <- TRUE
        if(check[1] != 0) beta[1:truth[1]] <- normalize(beta1) 
        if(check[2] !=0) beta[(truth[1]+1):length(beta_curr)] <- normalize(beta2) 
        break
      }
      beta <- c(normalize(beta1), normalize(beta2))
      conv <- (max(abs(beta - beta_curr)) < eps)
      beta_curr <- beta
    }
    geta <- z * B1 %*% delta1 + B2 %*% delta2
    g1 <- B1 %*% delta1
    g2 <- B2 %*% delta2
    loss <- -sum(log(1 + exp(geta))) + sum(y * geta)
    df1 <- sum(beta1!=0)
    df2 <- sum(beta2!=0)
    # df <- df1 + df2
    df <- df1 + df2 + 2*length(delta1)
    # bic <- -2 * loss + df * log(n) * log(log(p))
    bic <- -2 * loss + df * log(n) * log(p)
    aic <- -2 * loss + 2 * df
    # delta.var <- summary(fit.delta)$cov.unscaled  # not used 
    out <- list(beta1 = beta[1:truth[1]], beta2 = beta[(truth[1]+1):length(beta_curr)], B1 = B1, B2 = B2, delta1 = delta1, delta2 = delta2, iter = iter, g = geta, g1 = g1, g2 = g2, loss = loss, df = df, df1 = df1, df2 = df2, bic = bic, aic = aic, x = x, y = y, z = z, knots = knots, flag = flag, p=p, conv=conv, final.x = newx, final.off = off)
    class(out) <- "cste"
    return(out)
  }
}





