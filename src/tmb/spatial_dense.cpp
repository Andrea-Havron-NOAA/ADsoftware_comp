//spatial Poisson model
#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density; //namespace with GMRF function
  using namespace Eigen; //namespace with SparseMatrix declaration

  DATA_VECTOR(y);
  DATA_MATRIX(dd); //distance matrix
  DATA_VECTOR(prior_mean); 
  DATA_VECTOR(prior_sd); 
  DATA_INTEGER(prior_type);

  PARAMETER(b0); //intercept
  PARAMETER(ln_phi); //spatial correlation decay
  PARAMETER_VECTOR(omega); //spatial random effect

  int i;
  int n = y.size();

  Type nll = 0.0;

  if(prior_type == 1){
    nll -= dnorm(b0, Type(0), Type(5), true);
    nll -= dnorm(ln_kappa, prior_mean(0), prior_sd(0));
  }

  Type phi = exp(ln_phi);
  Type spsd = sqrt(exp(ln_spvar));
  //when spatial models are dense or nonseparable, it may be necessary to fix spatial parameters
  //Type marg_sp_sd = sqrt(0.75)
  //Type Range = 50;//exp(ln_phi);
  //Type kappa = sqrt(8)/Range;
  //Type tau =  1/(2*sqrt(M_PI)*kappa*marg_sp_sd);

  //Spatial Likelihood
  //Define precision matrix
  matrix<Type> cov(n,n); 
  for (i=0;i<n;i++)
  {
    cov(i,i)=Type(1);
    for ( j=0;j<i;j++)
    {
      cov(i,j)=exp(-phi*dd(i,j));			// Exponentially decaying correlation
      cov(j,i)=cov(i,j);
    }
  }
  
  MVNORM_t<Type> neg_log_density(cov);
  res+=neg_log_density(omega);

  //Data Likelihood

  vector<Type> eta(n);

  for(i=0; i<n; i++){
    eta(i) = b0 + omega(i);
    nll -= dpois(y(i), exp(eta(i)), true);
  }

  SIMULATE{
    y = rpois(exp(eta(i)));
    REPORT(y);
  }


  REPORT(omega);
  REPORT(nll);

  return nll;

}
