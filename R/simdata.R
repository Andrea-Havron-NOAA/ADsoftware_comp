## Simulate data FIMS architechture comparison
## Currently three models:
## 1. State Space Gompertz
## 2. State Space Logistic growth model
## 3. Spatial Poisson

library(mvtnorm)

gendat <- function(seed = NULL, N=100,
                   theta = c(2,.8), #gompertz: c(alpha=2,beta=.8); logistic:c(r=.2,K=100); spatial: c(b0 = 2,Range=50,SpVar=.75)
                   u1 = 4,
                   var = list(proc=0.1, obs=0.5), #gompertz: proc=0.1,obs=0.5; logistic: proc=0.01, obs=0.001
                   mod.name,
                   prior.type = 0,
                   mspp.dim = NULL,
                   spp.b0 = NULL,
                   loadings = NULL
                   ){
 
  if(mod.name == 'gompertz'){
    set.seed(seed)
    eta <- u <-  rep(0,N)
    u[1] <- u1
    for(i in 2:N){
      eta[i] <- theta[1] + theta[2]*u[i-1]
      u[i] <- rnorm(1, eta[i], sqrt(var$proc))
    }
    y <- rnorm(N,u,sqrt(var$obs))
  }
  
  if(mod.name == 'logistic'){

    r <- theta[1]
    K <- theta[2]
   
    set.seed(seed)
    eta <- u <-  rep(0,N)
    u[1] <- u1
    for(i in 2:N){
      eta[i] <- log(u[i-1] + r*u[i-1]*(1-u[i-1]/K))
      u[i] <- rlnorm(1, eta[i], sqrt(var$proc))
    }
    y <- rlnorm(N,log(u),sqrt(var$obs))
  }
  
  if(mod.name == 'spatial'){
    set.seed(seed)
    grid.xy <- expand.grid(x=seq(0,100,1), y=seq(0,100,1))
    d <- as.matrix(dist(grid.xy, upper = TRUE))
    u <- as.vector(sim.omega(theta[2], theta[3], d))
    y.sim <- rpois(nrow(grid.xy), exp(theta[1]+u))
    
   # samp.idx <- sample(1:nrow(grid.xy), N)
   # y <- data.frame(grid.xy[samp.idx,], z=y.sim[samp.idx])
    y <- data.frame(grid.xy, z=y.sim, u=u)
  }
  
  if(mod.name == 'spatial_mspp'){ 
    set.seed(seed)
    grid.xy <- expand.grid(x=seq(0,100,1), y=seq(0,100,1))
    d <- as.matrix(dist(grid.xy, upper = TRUE))
    nf <- mspp.dim[1]
    nj <- mspp.dim[2]
    u <- matrix(0, nrow(grid.xy), nf)
    y.sim <- matrix(0, nrow(grid.xy), nj)
    for(f in 1:nf){
      u[,f] <- as.vector(sim.omega(theta[1], theta[2], d))
    }
    L <- matrix(0, nj, nf)
    cnt <- 1
    for(f in 1:nf){
      for(j in 1:nj){
        if(j>f){
          L[j,f] <- loadings[cnt]
          cnt <- cnt + 1
        } else {
          L[j,f] <- 0
        }
      }
    }
    U <- u %*% t(L)
    for(j in 1:nj){
      y.sim[,j] <- rpois(nrow(grid.xy), exp(spp.b0[j]+U[,j]))
    }
    
    y <- data.frame(grid.xy, z=y.sim, U=U)
  }
  
  
  return(y)
}

## function to simulate spatial MVN using a Matern covariance function
sim.omega <- function(Range, sig2, Dmat, Nu = 1){
  Kappa <- sqrt(8)/Range
  N <- dim(Dmat)[1]
  Sigma <- matrix(NA,N,N)
  #matern covariance
  Sigma <- sig2 * 2^(1-Nu)/gamma(Nu) * (Kappa * Dmat) * besselK(Kappa * Dmat, Nu)
  diag(Sigma) <- sig2
  omega <- rmvnorm(1, rep(0,N), Sigma, method = "chol")
  
  return(omega)
}
