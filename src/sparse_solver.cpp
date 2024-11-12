#include <Rcpp.h>
using namespace Rcpp;

// calculate eta  
NumericVector etaC(NumericMatrix x, double a, NumericVector beta, NumericVector off){ 
  int n = x.nrow();
  int p = x.ncol();
  NumericVector ans(n);
  for (int i=0; i<n; ++i){
    for (int j=0; j<p; ++j){
      ans[i] += x(i, j)*beta[j];
    }
    ans[i] += off[i] + a;
  }
  return ans; 
}

// Soft thresholding operator
double st(double z, double lam) {
  if (z > lam)  return z - lam ;
  else if (fabs(z) <= lam) return 0;
  else return z + lam;
}

// The univariate solution for a MCP-penalized regression
double mcp(double z, double lam, double gam) {
  if (fabs(z) <= lam*gam) return st(z, lam)/(1-1/gam);
  else return z;
}

// The univariate solution for a SCAD-penalized regression
double scad(double z, double lam, double gam) {
  if (fabs(z) <= 2*lam) return st(z, lam);
  else if (fabs(z) > gam*lam) return z;
  else return st(z, gam*lam/(gam-1))/(1-1/(gam-1));
}

// check for convergence
double checkC(NumericVector b, NumericVector bold, double eps) {
  int converged = 1;
  int p = b.size();
  for (int i = 0; i < p; ++i) {
    if(fabs(b[i]-bold[i]) > eps) {
      converged = 0;
      break;
    }
  }
  return converged;
}

double cross(NumericMatrix x, NumericVector w, NumericVector r, int j, double ind) {
  double total = 0;
  int n = x.nrow();
  if(ind){
    for(int i=0; i<n; i++) total += w[i]*pow(x(i, j), 2)/n;
  } else {
    for(int i=0; i<n; i++) total += w[i]*x(i, j)*r[i]/n;
  }
  return total;
}

//' Solve the penalized logistic regression.
//' @param y samples of binary outcome which is a \eqn{n*1} vector.
//' @param x samples of covariates which is a \eqn{n*p} matrix.
//' @param off offset in logistic regression.
//' @param pen 1: MCP estimator; 2: SCAD estimator.
//' @param lam value of the lasso penalty parameter \eqn{\lambda} for \eqn{\beta_1} and \eqn{\beta_2}.
//' @param beta initial estimates.
//' @return A numeric vector, estimate of beta
//' @export 
// [[Rcpp::export]]
NumericVector penC(NumericMatrix x, NumericVector y, NumericVector off, NumericVector beta, double lam, double pen) {
  int n = x.nrow();
  int p = x.ncol();
  double pi = 0;
  double iter = 0;
  double shift = 0;
  double eps = 0.001;
  NumericMatrix xx = clone(x);
  NumericVector beta_est(p);
  NumericVector z(p);
  NumericVector eta(n);
  NumericVector w(n);
  NumericVector r(n);
  NumericVector v(p);
  NumericVector betaold(p);
  NumericVector avg(p);
  NumericVector std(p);
  
  // column mean, column standard deviation, normalize data 
  for(int j=0; j<p; j++) {
    avg[j] = sum( x(_, j) ) / n;
    std[j] = sqrt( sum( (x(_, j) - avg[j]) * (x(_, j) - avg[j]) ) / n );
    xx(_, j) = (x(_, j) - avg[j]) / std[j];
  }
  
  //initial intercept
  double a = 0;
  a = sum(avg * beta / std);
  
  while(iter < 1000) {
    // Previous iteration
    for (int j=0; j<p; j++) {
      betaold[j] = beta[j];
    }
    
    eta = etaC(xx, a, beta, off);
    a = 0;          // why? 
    
    // Calculate w, r
    for (int i=0; i<n; i++) {
      if (eta[i] > 10) {
        pi = 1;
        w[i] = .0001;
      } else if (eta[i] < -10) {
        pi = 0;
        w[i] = .0001;
      } else {
        pi = exp(eta[i])/(1+exp(eta[i]));
        w[i] = pi*(1-pi);
      }
      r[i] = (y[i] - pi)/w[i];
      
    }
    // Covariates
    for (int j=0; j<p; j++) {
      
      // Calculate v_j
      // v[j] = cross(xx, w, r, j, 1);
      v[j] =  sum(w * xx(_, j)  * xx(_, j)) / n;  
      
      // Update z_j
      // z[j] = cross(xx, w, r, j, 0)+v[j]*beta[j];
      z[j] =  sum(w * xx(_,j) * r) / n +v[j]*beta[j]; 
      
      
      // Update b_j
      if (pen == 1){
        beta[j] = mcp(z[j], lam, 3.0)/v[j];
      } else{
        beta[j] = scad(z[j], lam, 3.7)/v[j];
      }
      
      // Update r
      shift = beta[j]-betaold[j];
      for (int i=0; i<n; i++) {
        r[i] -= shift*xx(i,j);
      }
    }
    // Intercept
    a = sum(avg * beta / std);
    
    // Check for convergence
    if(checkC(beta, betaold, eps)){
      break;
    }
    iter += 1;
  }
  for (int j=0; j<p; j++) {
    beta_est[j] = beta[j]/std[j];
  }
  return beta_est;
}
