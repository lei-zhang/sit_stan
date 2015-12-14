source("_scripts/stan_source.R")
library(ggplot2)
library(reshape2)

fit_RL    <- readRDS("_outputs/fit_RL.RData")
lik_RL    <- extract_log_lik(fit_RL, parameter_name = "log_lik")
loo_RL    <- loo(lik_RL, cores = 4)
loo_RL_pw <- as.matrix(loo_RL$pointwise[,3] )

fit_RLnc    <- readRDS("_outputs/fit_RLnc_2lr.RData")
lik_RLnc    <- extract_log_lik(fit_RLnc, parameter_name = "log_lik")
loo_RLnc    <- loo(lik_RLnc, cores = 4)
loo_RLnc_pw <- as.matrix(loo_RLnc$pointwise[,3] )

fit_RLcoh    <- readRDS("_outputs/fit_RLcoh_2lr_cfa.RData")
lik_RLcoh    <- extract_log_lik(fit_RLcoh, parameter_name = "log_lik")
loo_RLcoh    <- loo(lik_RLcoh, cores = 3)
loo_RLcoh_pw <- as.matrix(loo_RLcoh$pointwise[,3] )

fit_RLcr    <- readRDS("_outputs/fit_RLcumrew_2lr.RData")
lik_RLcr    <- extract_log_lik(fit_RLcr, parameter_name = "log_lik")
loo_RLcr    <- loo(lik_RLcr, cores = 3)
loo_RLcr_pw <- as.matrix(loo_RLcr$pointwise[,3] )

save(loo_RL_pw, loo_RLnc_pw, loo_RLcoh_pw, loo_RLcr_pw,
     file = '_outputs/looic.RData')

looMat = cbind(loo_RL_pw, loo_RLnc_pw, loo_RLcoh_pw, loo_RLcr_pw)
save(looMat, file = '_outputs/looicMat.RData')

#### plot ####
load(file = '_outputs/looicMat.RData')

looDF <- melt(looMat)
colnames(looDF) <- c('subID', 'modelName', 'LOOIC')
looDF$modelName <- c(rep('RL', 129), rep('RLnc',129), 
                     rep('RLcoh',129), rep('RLcumrew',129))
modelName <- c('RL','RLnc','RLcoh','RLcumrew')
looDF$modelName <- factor(looDF$modelName, levels = modelName)
theme_set(theme_gray(base_size = 20))

ggplot(data=looDF, aes(modelName, LOOIC)) + geom_boxplot(aes(fill = modelName)) 

#### individual ####
source("http://bioconductor.org/biocLite.R")
biocLite("Biobase")

minlooic <- do.call(pmin, as.data.frame(looMat))
minMat <- kronecker(matrix(1,1,4),minlooic)
colSums((minMat == looMat))

#### use subj66 for plot the ppc
source('_scripts/ppc.R')
ppc_RL    <- ppc(fit_RL, sid = 66)
ppc_RLnc  <- ppc(fit_RLnc, sid = 66)
ppc_RLcoh <- ppc(fit_RLcoh, sid = 66)
ppc_RLcr  <- ppc(fit_RLcr, sid = 66)

save(ppc_RL,ppc_RLnc,ppc_RLcoh,ppc_RLcr,file= 'ppc2.RData')













