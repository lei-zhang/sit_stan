beta_barplot <- function(stanfit, six = FALSE) {
# beta_barplot() plots the BETAs' posterior group mean with a 95% HDI
  
  library(ggplot2)
  library(coda)
  
  if ( class(stanfit)=='stanfit' ) {
    stanfit <- stanfit
  } else {
    stanfit <- stanfit$fit
  }
  
  #### get BETAs from Stanfit output -------------------------------------------------------------
  parm   <- get_posterior_mean(stanfit, 'beta_mu')
  b_mean <- as.matrix(parm[,5])
  nb     <- length(b_mean)
  b_name <- paste('beta', 1:nb)
  
  df <- data.frame(b_mean = b_mean, b_name = b_name )

  #### extract mcmc and calculate the 95% HDI -----------------------------------------------------
  mcmcCoda = mcmc.list( lapply( 1:ncol(stanfit) , function(x) {mcmc(as.array(stanfit)[,x,])}))
  
  HDI <- array(0, c(nb,2))
  for (B in 1:nb) {
    HDI[B,] <- HDIofMCMC( mcmcCoda[, paste0("beta_mu[",B,"]") ] )
  }
  HDI_lb <- HDI[,1];    HDI_ub <- HDI[,2]
  df$HDI_ub <- HDI_ub;  df$HDI_lb <- HDI_lb
  df <- df[c('b_name','b_mean','HDI_ub','HDI_lb')]

  #### plot #### ---------------------------------------------------------------------------------
  theme_set(theme_bw(base_size = 18))
  
  g <- ggplot(df, aes(x = b_name, y = b_mean))
  dodge  <- position_dodge(width=0.9)
  limits <- aes(ymax = HDI_ub, ymin = HDI_lb)
  
  if ( six == FALSE) {
    g <- g + geom_bar(position=dodge, stat="identity", width = .8,
                      fill = "deepskyblue2", colour = "black")
    g <- g + geom_errorbar(limits, position=dodge, width=0.25, 
                           size = 1, colour = "deepskyblue4")
  } else {
    g <- g + geom_bar(position=dodge, stat="identity", width = .8,  
                      fill = c("#D55E00","#56B4E9","#999999","#999999","#009E73","#009E73"))
    g <- g + geom_errorbar(limits, position=dodge, width=0.25, size = 1,
                           colour = c("#D55E00","#56B4E9","#999999","#999999","#009E73","#009E73"))
  }
                         
  g <- g + geom_hline(yintercept=0)
  g <- g + ylab("beta values") + xlab("")
  print(g)
  
  return(df)
}


#### HDI function adapted from Kruschke's book (2015)

HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  
  if ( class(sampleVec) == "mcmc.list" ) {
    sampleVec = as.matrix(sampleVec)
  }
  
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}
