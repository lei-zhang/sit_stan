beta_corplot <- function( stanfit, splt = NULL ) {
# beta_corplot() plots the correlation matrix of individuals' BETAs
  
  library(corrplot)
  
  if ( class(stanfit)=='stanfit' ) {
    stanfit <- stanfit
  } else {
    stanfit <- stanfit$fit
  }
  
  #### get BETAs from Stanfit output -------------------------------------------------------------
  ns <- 129
  parm <- get_posterior_mean(stanfit, 'beta')
  bMat <- matrix(parm[,5], nrow = ns)
  nb   <- dim(bMat)[2]
  colnames(bMat) <- paste('beta', 1:nb)
  
  ### get all the parameters ----------------------
#   parm <- get_posterior_mean(stanfit, c('lr', 'disc', 'beta'))
#   parMat <- matrix(parm[,5], nrow = ns)  
#   colnames(parMat) <- c('lr', 'disc',paste('beta', 1:6))
#   writeMat('_outputs/param_beta4_alt6.mat', parMat = parMat)


  ------------------------------------------------------------------------  
  if (is.null(splt) ) {
    corrMat <- cor(bMat)
    corrplot(corrMat, method = 'square', diag = FALSE, addCoef.col = T)
  } else if (splt == 6) {
    bMat_pos <- bMat[bMat[,6] >= 0,]
    bMat_neg <- bMat[bMat[,6] < 0,]
    corrMat_pos <- cor(bMat_pos)
    corrMat_neg <- cor(bMat_neg)
    par(mfrow = c(1,2))
    corrplot(corrMat_pos, method = 'square', diag = FALSE, addCoef.col = T)
    corrplot(corrMat_neg, method = 'square', diag = FALSE, addCoef.col = T)
  } 
}