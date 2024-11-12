#' Calculate simultaneous confidence bands of CSTE curve for binary outcome.  
#' 
#' This function calculates simultaneous confidence bands of CSTE curve for binary 
#' outcome. 
#' 
#' 
#'@param x samples of predictor, which is a \eqn{m*p} matrix.
#'@param fit a S3 class of cste.
#'@param h kernel bandwidth.
#'@param alpha the simultaneous confidence bands are of \eqn{1-\alpha} confidence level.
#'
#'@return A list which includes:
#' \itemize{
#'    \item \code{or_x}: the ordered value of \eqn{X\beta_1}. 
#'    \item \code{fit_x}: the fitted value of CSTE curve corresponding to \code{or_x}.
#'    \item \code{lower_bound}: the lower bound of CSTE's simultaneous confidence band.
#'    \item \code{upper_bound}: the upper bound of CSTE's simultaneous confidence band.
#' }
#'
#'
#'
#' @references
#' Guo W., Zhou X. and Ma S. (2021).
#' Estimation of Optimal Individualized Treatment Rules
#' Using a Covariate-Specific Treatment Effect Curve with 
#' High-dimensional Covariates,
#' \emph{Journal of the American Statistical Association}, 116(533), 309-321
#' 
#' @seealso \code{\link{cste_bin}}

cste_bin_SCB <- function(x, fit, h = NULL, alpha = 0.05){
  u1 <- pu(x, fit$beta1)$u
  u2 <- pu(x, fit$beta2)$u
  sbk <- prev_fit_cste(u1, u2, fit)
  if(is.null(h)){
    # optimal bandwidth 
    # h <- sbk$h * (log(sbk$n))^(-0.25)
    h <- sbk$h    
  }
  newx <- seq(min(u1)+h, max(u1)-h, length = 100)
  fit.x <- sapply(newx, function(xx) coef(glm(sbk$y~1, weights=dnorm((xx - sbk$u1)/h)/h ,
                        family=quasibinomial(link = "logit"), offset=sbk$fit_g2)))
  
  # estimate sigma_b^2(x)
  fit_sb <- predict(sbk$fit_sigma_b, newx)$y
  
  # estimate sigma^2(x)
  fit_s <- predict(sbk$fit_sigma_x, newx)$y
  
  # calculate inflation factor 
  # alpha <- 0.05
  Ck_d <- 1/(2*sqrt(pi))
  Ck_n <- 1/(4*sqrt(pi))
  Ck <- Ck_n/Ck_d
  mu2k <- 1/(2*sqrt(2))
  ah <- sqrt(-2 * log(h))
  Qh <- ah + (log(sqrt(Ck)/(2*pi)) - log(-log(sqrt(1 - alpha))))/ah
  
  # estimate density of u1
  h_x <- bw.nrd0(sbk$u1)
  f_x <- sapply(newx, function(xx) mean(dnorm((xx - u1)/h_x)/h_x))
  
  # calculate variance
  D <- fit_sb * f_x
  v_sq <- Ck_d * f_x * fit_s
  id_rm <- D < 0 | v_sq < 0
  
  g_sigma <- sbk$n^(-0.5) * h^(-0.5) * sqrt(v_sq[!id_rm]) / D[!id_rm]
  
  # calculate SCC
  L <- fit.x[!id_rm] - Qh * g_sigma
  U <- fit.x[!id_rm] + Qh * g_sigma
  
  # without model selection h , 0.006
  or_x <- pu_inv(x, fit$beta1, newx)
  # mycol <- gg_color(2)
  # or_x <- quantile(seq(min(x%*%fit$beta1), max(x%*%fit$beta1),length=1000), newx)
  
  return(list(or_x = or_x[!id_rm],fit_x = fit.x[!id_rm],
              lower_bound = L, upper_bound = U))
}