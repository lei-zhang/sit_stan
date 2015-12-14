get_mcmc <- function(stanFit){
  
  library(coda)     # for mcmc.list()
  library(runjags)  # for combine.mcmc()
  
  mcmcCoda <- mcmc.list( lapply( 1:ncol(stanFit), function(x){mcmc(as.array(stanFit)[,x,])}))
  mcmc     <- combine.mcmc(mcmcCoda)
  
  return(mcmc)
}
