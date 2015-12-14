source("_scripts/stan_run.R")
source("_scripts/compute_map.r")
source("_scripts/get_mcmc.R")
source("_scripts/cal_prob_v1.R")
source("_scripts/cal_prob_v2.R")
source("_scripts/ppc.R")
source("_scripts/corr_b1b2.R")
source("_scripts/extract_looic.R")
source("_scripts/extract_otheric.R")
source("_scripts/beta_barplot.R")
source("_scripts/beta_corplot.R")
#source("_scripts/DBDA2E-utilities.R")

library(rstan)
library(shinystan)
library(loo)
library(R.matlab)


# ---- convert .mat data into .RData
# t <- readMat("_data/data3_129.mat")
# mydata <- t$data3
# save(mydata, file = '_data/sit_reversal_betnoshow_129.rdata')

# ---- obtain mcmc -----------------------
# mcmcCoda <- mcmc.list( lapply( 1:ncol(out1), function(x){mcmc(as.array(out1)[,x,])}))
# mcmc     <- combine.mcmc(mcmcCoda)
# 
# # ---- diagnostic plots ------------------
# diagMCMC(mcmcCoda, parName=c("lr_mu"))
# diagMCMC(mcmcCoda, parName=c("tau_mu"))
# diagMCMC(mcmcCoda, parName=c("lr[1]"))
# diagMCMC(mcmcCoda, parName=c("tau[1]"))
# 
# plotPost(mcmcCoda[,"lr_mu"], main = "lr_mu", xlab = bquote(lr_mu), cenTend = "mode")
# plotPost(mcmcCoda[,"tau_mu"], main = "tau_mu", xlab = bquote(tau_mu), cenTend = "mode")
# plotPost(mcmcCoda[,"lr[1]"], main = "lr[1]", xlab = bquote(lr[1]), cenTend = "mode")
# plotPost(mcmcCoda[,"tau[1]"], main = "tau[1]", xlab = bquote(tau[1]), cenTend = "mode")

## end of script ##