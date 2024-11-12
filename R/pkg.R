#' @useDynLib CSTE
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @importFrom fda bsplineS 
#' @import Rcpp
#' @importFrom splines bs
#' @importFrom stats optim approx binomial coef dbeta density glm.fit pbeta predict qbeta quantile smooth.spline bw.nrd0 complete.cases dnorm glm rnorm
#' @importFrom stats na.omit qnorm quasibinomial
#' @importFrom survival coxph Surv survfit basehaz 
#' @importFrom locpol gaussK
#' @importFrom dfoptim nmk


NULL