corr_b1b2 <- function(stanfit, plot_sct=FALSE) {
  
  if ( class(stanfit)=='stanfit' ) {
    stanfit <- stanfit
  } else {
    stanfit <- stanfit$fit
  }
  
  parm <- get_posterior_mean(stanfit, 'beta')
  parm <- as.matrix(parm[1:258,5])
  b1 <- parm[1:129]
  b2 <- parm[130:258]
  
  r <- cor.test(b1,b2)
  
  if (plot_sct == TRUE) {
    library(ggplot2)
    q <-  qplot(b1,b2,geom='point')
    print(q)
  }
    
  return(r)
}
  