setupADMB <- function(path, modName, doRE){ 
  wd <- getwd()
  setwd(path)
  compile_admb(modName, re = doRE, verbose = TRUE)
  setwd(wd)
}