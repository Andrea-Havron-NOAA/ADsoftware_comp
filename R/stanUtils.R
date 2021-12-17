#' Make STAN data list
#'
#' @param obs vector of observations
#' @param ptype integer specifying if using improper priors (ptype=0) 
#'              or vague proper priors (ptype=1)
#'
#' @return data list
#' @export
mkSTANLogisticDat <- function(obs,
                      hyperParms = list(
                        hyperSig = c(0.1),
                        hyperTau = c(0.1),
                        hyperTheta1 = c(-1,4),
                        hyperTheta2 = c(5,4)
                      )){
  return(list(N = length(obs), 
              y = obs,
              hyp_sig = hyperParms$hyperSig,
              hyp_tau = hyperParms$hyperTau,
              hyp_theta1 = hyperParms$hyperTheta1,
              hyp_theta2 = hyperParms$hyperTheta2,
              prior_type = 0))
  
}


mkSTANSpatialInits <- function(df,pr,method){
  Loc <- df[,1:2]
  y <- df[,3]
  
  #build INLA mesh
  mesh <- inla.mesh.2d(Loc, max.edge = c(10,20), offset = c(5,25))
  #calculate sparse distance matrix components
  spde <- inla.spde2.matern(mesh)
  
  Dat <- list(N = length(y),
              NV = mesh$n,
              y = y,
              vi = mesh$idx$loc,
              M0 = as.matrix(spde$param.inla$M0),
              M1 = as.matrix(spde$param.inla$M1),
              M2 = as.matrix(spde$param.inla$M2),
              #lnkapPr = c(0,0),
              #lntauPr = c(0,0),
              kap_tau_pr_mu = c(0,0),
              kap_tau_pr_var = c(0,0),
              prior_type = 0,
              kappa = sqrt(8)/50,
              tau = 1/(sqrt(4*pi*0.75)*sqrt(8)/50))
  if(pr == 1){
    #Dat$lnkapPr <- c(max(dist(Loc))*.2, 10)
    #Dat$lntauPr <- c(0,2)
    #Dat$kap_tau_pr_mu <- c(log(max(dist(Loc))*.2), 0)
    #Dat$kap_tau_pr_var <- c(10, 2)
    Dat$prior_type = 1
  }

  init.list <- list(Dat = Dat)
  return(init.list)
}


