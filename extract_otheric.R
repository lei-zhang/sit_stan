extract_otheric <- function(stanfit, choice=2, np) {
  # extract AIC, NIC, WAIC, from mcmc object
  # np: number of parameters
  
  L <- list()
  
  if ( class(stanfit)=='stanfit' ) {
    stanfit <- stanfit
  } else {
    stanfit <- stanfit$fit
  }
  
  if ( is.null(choice) )  {
    par_name <- "log_lik"
  } else if (choice == 1) {
    par_name <- "log_likc1"
  } else if (choice == 2) {
    par_name <- "log_likc2"
  }
  
  lik     <- extract_log_lik(stanfit, parameter_name = par_name)
  myWAIC  <- waic(lik)
  
  subjLik <- colMeans(lik) # estimated loglikelihood per subject from the posterior loglik matrix
  ns <- dim(lik)[2] # number of subjects
  myAIC <- -2 * subjLik + 2 * np
  myBIC <- -2 * subjLik + log(ns) * np
  
  L$myAIC  <- sum(myAIC)
  L$myBIC  <- sum(myBIC)
  L$myWAIC <- myWAIC
  
  return(L)
}