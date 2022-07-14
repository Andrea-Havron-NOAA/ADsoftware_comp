#' Setup C++ TMB code 
#'
#' @param comp boolean that determines whether or not to compile the C++ file
#'
#' @return DLL is compiled and loaded
#' @export
setupTMB <- function(dll.name, comp=FALSE){
  if(!(paste0(dll.name, '.dll') %in% list.files('src/tmb')) |
     comp==TRUE){
    try(dyn.unload(dynlib(paste0('src/tmb/', dll.name))))
    TMB::compile(paste0('src/tmb/', dll.name,'.cpp'))
  }
  suppressWarnings(dyn.load(dynlib(paste0('src/tmb/', dll.name))))
}

#'
#' @param obs vector of observation data
#' @param modname string specifying 'gompertz' or 'logistic'
#' @param prType Boolean indicating whether or not to include priors
#'
#' @return data list
#' @export
#'
mkTMBLogisticDat <- function(obs, modname, prType){
  dat <- list(y=obs)
  if(prType == 0){
    dat$hyperpars <- 0
  } else {
    #TODO: hyperpars need to be fixed in TMB cpp files
    dat$hyperpars <- c(-1,4,5,4,0.1,0.1) 
  }
  return(dat)
}

mkTMBSpatialInits <- function(df,pr,method, nf = NULL){
  Loc <- df[,1:2]
  y <- df[,3:ncol(df)]
  
  #build INLA mesh
  mesh <- inla.mesh.2d(Loc, max.edge = c(10,20), offset = c(5,25))
  #calculate sparse distance matrix components
  spde <- inla.spde2.matern(mesh)

  Dat <- list(y = y,
              v_i = mesh$idx$loc-1,
              M0 = spde$param.inla$M0,
              M1 = spde$param.inla$M1,
              M2 = spde$param.inla$M2,
              #kap_tau_pr_mu = c(0,0),
              #kap_tau_pr_var = c(0,0),
              prior_type = 0,
              kappa = sqrt(8)/50,
              tau = 1/(sqrt(4*pi*0.75)*sqrt(8)/50))
  #hyperpars = matrix(0,2,2))
  if(pr == 1){
    # Dat$hyperpars[1,] <- c(max(dist(Loc))*.2, 10)
    # Dat$hyperpars[2,] <- c(0,2)
    # Dat$kap_tau_pr_mu <- c(log(max(dist(Loc))*.2), 0)
    # Dat$kap_tau_pr_var <- c(10, 2)
    Dat$prior_type = 1
  }
  Par.fn <- function(){
    list(b0 = 0,#ln_phi = 0, ln_spvar = 0,
         omega = rep(0,mesh$n))
  }
  
  if(method == 'tmb_mspp'){
    Dat <- list(y = as.vector(as.matrix(y)), 
                v_i = rep(mesh$idx$loc-1, ncol(y)), 
                M0 = spde$param.inla$M0,
                M1 = spde$param.inla$M1,
                M2 = spde$param.inla$M2,
                ni = nrow(y),
                nj = ncol(y),
                nf = nf,
                prior_type = 0)
    Dat$prior_type <- pr
    nl <- ncol(y)*nf-(nf*(nf-1))/2
    Par.fn <- function(){
      list(b0 = rep(0,ncol(y)),
           loadings = rep(1,nl),
           omega = rep(0,mesh$n*nf))
    }
  } 
 
  
  init.list <- list(Dat = Dat, Par = Par.fn)
  return(init.list)
}

#' Helper function to run TMB code and return output

#' @param obj.args List of arguments for TMB MakeADFun() function
#' @param opt.args List of arguments for nlminb() function
#' @param control List controlling model runs and standard error reporting
#'
#' @return Fitted objective function, nlminb output, reported values from model, sdreport if true
#'
fitTMB <- function(obj.args, opt.args = list(control = list(iter = 800, eval = 800),
                                              scale = 1,
                                              lower = -Inf, upper = Inf ),
                    control = list(run.model = TRUE, do.sdreport = TRUE)){
  obj <- do.call(MakeADFun, obj.args)
  if(control$run.model){
    opt <- with(obj, do.call(nlminb,  c(list(par, fn, gr), opt.args) ))
    report <- obj$report(obj$env$last.par.best)
    if(control$do.sdreport){
      sdr <- sdreport(obj)
      fit.results <- list(obj = obj, opt = opt, report = report, sdr = sdr, aic = aic)
    } else {
      fit.results <- list(obj = obj, opt = opt, report = report, aic = aic)
    }
  } else {
    fit.results <- list(obj = obj, report = obj$report())
  }
  return(fit.results)
}