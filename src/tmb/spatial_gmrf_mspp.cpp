//spatial Poisson model
#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density; //namespace with GMRF function
  using namespace Eigen; //namespace with SparseMatrix declaration

  DATA_VECTOR(y);
  DATA_IVECTOR(v_i); //location index
  DATA_SPARSE_MATRIX(M0); //sparse distance matrices
  DATA_SPARSE_MATRIX(M1); //sparse distance matrices
  DATA_SPARSE_MATRIX(M2); //sparse distance matrices
  DATA_INTEGER(ni); //number of sites
  DATA_INTEGER(nj); //number of species
  DATA_INTEGER(nf); //number of factors
  DATA_INTEGER(prior_type);

  PARAMETER_VECTOR(b0); //species-specific intercept
  PARAMETER_VECTOR(loadings);
  PARAMETER_VECTOR(omega); //spatial random effect

  int n = y.size();
  int nv = M0.col(0).size();

  Type nll = 0.0;
  vector<Type> resid(2);

  if(prior_type == 1){
    nll -= sum(dnorm(b0, Type(0), Type(5), true));
    nll -= sum(dnorm(loadings, Type(0), Type(2), true));
  }

  Type marg_sp_sd = sqrt(0.75);//sqrt(exp(ln_spvar));
  Type Range = 50;//exp(ln_phi);
  Type kappa = sqrt(8)/Range;
  Type tau =  1/(2*sqrt(M_PI)*kappa*marg_sp_sd);

  //Spatial Likelihood
  //Define precision matrix
  int cnt = 0;
  matrix<Type> L(nj, nf);
  matrix<Type> Omega_vf(nv,nf);
  SparseMatrix<Type> Q = pow(kappa,4) * M0 + 2 * pow(kappa,2) * M1 + M2;
  for(int f=0; f<nf; f++){
    for(int j=0; j<nj; j++){
      if(j>=f){
        L(j,f) = loadings(cnt);
        cnt ++;
      } else {
        L(j,f) = 0.0;
      }
    }
    for(int v=0; v<nv; v++){
      int idx = f*nv+v;
      Omega_vf(v,f) = omega(idx);
    }
    nll += SCALE(GMRF(Q), 1/tau)(Omega_vf.col(f));
  }
  
  
  SIMULATE{
    for(int f=0; f<nf; f++){
      Omega_vf.col(f) = GMRF(Q).simulate()/tau;
    }
    REPORT(Omega_vf);
  }

  //Data Likelihood
  matrix<Type> Omega_vj = Omega_vf * L.transpose();
  vector<Type> eta(n);

  for(int i=0; i<ni; i++){
    for(int j=0; j<nj; j++){
      int idx = i*nj+j;
      eta(idx) = b0(j) + Omega_vj(v_i(idx),j);
    }
  }
  
  nll -= sum(dpois(y, exp(eta), true));

  SIMULATE{
    y = rpois(exp(eta));
    REPORT(y);
  }

  REPORT(Range);
  ADREPORT(Range);
  REPORT(marg_sp_sd);
  ADREPORT(marg_sp_sd);
  REPORT(omega);
  REPORT(nll);

  return nll;

}
