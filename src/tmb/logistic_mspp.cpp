#include <TMB.hpp>

// dlnorm
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator()()
{
  DATA_MATRIX(y);
  DATA_MATRIX(X);
  DATA_VECTOR(hyperpars); //theta1 mu, theta1 sd, theta2 mu, theta2 sd, sigma, tau

  PARAMETER_MATRIX(beta)
  PARAMETER_VECTOR(ln_K); vector<Type> K = exp(ln_K);
  PARAMETER_VECTOR(ln_sig);  vector<Type> sigma = exp(ln_sig);
  PARAMETER_VECTOR(ln_tau);  vector<Type> tau = exp(ln_tau);

  PARAMETER_MATRIX(u);

  int nt = y.col(0).size();
  int nj = y.row(0).size();
  int t;
  int j;
  Type nll = 0;


  matrix<Type> eta(nt,nj);
  if(hyperpars.size() > 1){
    nll -= sum(dlnorm(beta, hyperpars(0), hyperpars(1), true));
    nll -= sum(dlnorm(K, hyperpars(2), hyperpars(3), true));
    nll -= sum(dexp(sigma, hyperpars(4), true));
    nll -= sum(dexp(tau, hyperpars(5), true));
  }

  matrix<Type> r = exp(X * beta); //would invlogit be better?
  for (t = 1; t < nt; t++) {
    for (j = 0; j < nj; j++) {
        eta(t,j) = log(u((t - 1), j) + r(t,j) * u((t - 1), j) * (1 - u((t - 1), j) / K(j)));
        nll -= dlnorm(u(t, j), eta(t,j), sigma(j), true);
    }
  }
  for (int t = 0; t < nt; t++) {
    for (int j = 0; j < nj; j++) {
        nll -= dlnorm(y(t,j), log(u(t, j)), tau(j), true);
    }
  }
    
  REPORT(r);
  REPORT(K);
  ADREPORT(r);
  ADREPORT(K);
  REPORT(sigma);
  REPORT(tau);
  ADREPORT(sigma);
  ADREPORT(tau);
  
  return(nll);
}
