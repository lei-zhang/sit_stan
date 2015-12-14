extract_looic <- function(stanfit, choice=2, core=4) {
  # extract LOOIC, from mcmc object
  
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
  
  lik    <- extract_log_lik(stanfit, parameter_name = par_name)
  looic  <- loo(lik, cores = core)
  
  return(looic)
}