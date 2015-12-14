compute_map <- function(fitObj, meth = "density", calBeta = F) {
  
  library(modeest); library(coda); library(runjags)
  
  startTime = Sys.time()
  
  if ( class(fitObj) == 'mcmc' ) {
    mcmc <- fitObj
  } else if ( class(fitObj) == 'mcmc.list' ) {
    mcmc <- combine.mcmc(fitObj)
  } else {
    if (class(fitObj) != 'stanfit') {
      fitObj <- fitObj$fit
    }
    mcmcCoda <- mcmc.list( lapply( 1:ncol(fitObj) , function(x) {mcmc(as.array(fitObj)[,x,])}))
    mcmc     <- combine.mcmc(mcmcCoda)
  }
  
  if (calBeta == T) {
    b_start <- which( colnames(mcmc) == 'beta[1,1]' )
    b_end   <- which( colnames(mcmc) == 'beta[6,129]' )
    mcmc    <- mcmc[, b_start:b_end]
    sz <- dim(mcmc) 
    map <- array(0,dim=c(sz[2],1))
  } else {
    sz <- dim(mcmc) 
    map <- array(0,dim=c(sz[2],1))
  }
  
  for (f in 1:sz[2]) {
    tmp <- mlv(mcmc[,f], method=meth)
    map[f] <- mean(tmp$M)
  }
  rownames(map) <- colnames(mcmc)
  
  if (calBeta == T) {
    map <- t(matrix(map, nrow = 6))
  }
  
  endTime = Sys.time()
  cat("It took",as.character.Date(endTime - startTime), "\n")
  return(map)
}
