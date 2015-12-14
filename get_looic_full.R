source("_scripts/stan_source.R")
library(ggplot2)
library(reshape2)

fit_RL    <- readRDS("_outputs/fit_RL.RData")
lik_RL    <- extract_log_lik(fit_RL, parameter_name = "log_lik")
loo_RL    <- loo(lik_RL, cores = 4)
loo_RL_pw <- as.matrix(loo_RL$pointwise[,3] )

fit_RLnc    <- readRDS("_outputs/fit_RLnc.RData")
lik_RLnc    <- extract_log_lik(fit_RLnc, parameter_name = "log_lik")
loo_RLnc    <- loo(lik_RLnc, cores = 4)
loo_RLnc_pw <- as.matrix(loo_RLnc$pointwise[,3] )

fit_RLnc_2lr    <- readRDS("_outputs/fit_RLnc_2lr.RData")
lik_RLnc_2lr    <- extract_log_lik(fit_RLnc_2lr, parameter_name = "log_lik")
loo_RLnc_2lr    <- loo(lik_RLnc_2lr, cores = 4)
loo_RLnc_2lr_pw <- as.matrix(loo_RLnc_2lr$pointwise[,3] )

fit_RLnc_cfa    <- readRDS("_outputs/fit_RLnc_cfa.RData")
lik_RLnc_cfa    <- extract_log_lik(fit_RLnc_cfa, parameter_name = "log_lik")
loo_RLnc_cfa    <- loo(lik_RLnc_cfa, cores = 4)
loo_RLnc_cfa_pw <- as.matrix(loo_RLnc_cfa$pointwise[,3] )

fit_RLnc_2lr_cfa    <- readRDS("_outputs/fit_RLnc_2lr_cfa.RData")
lik_RLnc_2lr_cfa    <- extract_log_lik(fit_RLnc_2lr_cfa, parameter_name = "log_lik")
loo_RLnc_2lr_cfa    <- loo(lik_RLnc_2lr_cfa, cores = 4)
loo_RLnc_2lr_cfa_pw <- as.matrix(loo_RLnc_2lr_cfa$pointwise[,3] )

fit_RLcoh    <- readRDS("_outputs/fit_RLcoh.RData")
lik_RLcoh    <- extract_log_lik(fit_RLcoh, parameter_name = "log_lik")
loo_RLcoh    <- loo(lik_RLcoh, cores = 3)
loo_RLcoh_pw <- as.matrix(loo_RLcoh$pointwise[,3] )

fit_RLcoh_2lr    <- readRDS("_outputs/fit_RLcoh_2lr.RData")
lik_RLcoh_2lr    <- extract_log_lik(fit_RLcoh_2lr, parameter_name = "log_lik")
loo_RLcoh_2lr    <- loo(lik_RLcoh_2lr, cores = 3)
loo_RLcoh_2lr_pw <- as.matrix(loo_RLcoh_2lr$pointwise[,3] )

fit_RLcoh_cfa    <- readRDS("_outputs/fit_RLcoh_cfa.RData")
lik_RLcoh_cfa    <- extract_log_lik(fit_RLcoh_cfa, parameter_name = "log_lik")
loo_RLcoh_cfa    <- loo(lik_RLcoh_cfa, cores = 3)
loo_RLcoh_cfa_pw <- as.matrix(loo_RLcoh_cfa$pointwise[,3] )

fit_RLcoh_2lr_cfa    <- readRDS("_outputs/fit_RLcoh_2lr_cfa.RData")
lik_RLcoh_2lr_cfa    <- extract_log_lik(fit_RLcoh_2lr_cfa, parameter_name = "log_lik")
loo_RLcoh_2lr_cfa    <- loo(lik_RLcoh_2lr_cfa, cores = 3)
loo_RLcoh_2lr_cfa_pw <- as.matrix(loo_RLcoh_2lr_cfa$pointwise[,3] )

fit_RLcr    <- readRDS("_outputs/fit_RLcumrew.RData")
lik_RLcr    <- extract_log_lik(fit_RLcr, parameter_name = "log_lik")
loo_RLcr    <- loo(lik_RLcr, cores = 3)
loo_RLcr_pw <- as.matrix(loo_RLcr$pointwise[,3] )

fit_RLcr_2lr    <- readRDS("_outputs/fit_RLcumrew_2lr.RData")
lik_RLcr_2lr    <- extract_log_lik(fit_RLcr_2lr, parameter_name = "log_lik")
loo_RLcr_2lr    <- loo(lik_RLcr_2lr, cores = 3)
loo_RLcr_2lr_pw <- as.matrix(loo_RLcr_2lr$pointwise[,3] )

fit_RLcr_cfa    <- readRDS("_outputs/fit_RLcumrew_cfa.RData")
lik_RLcr_cfa    <- extract_log_lik(fit_RLcr_cfa, parameter_name = "log_lik")
loo_RLcr_cfa    <- loo(lik_RLcr_cfa, cores = 3)
loo_RLcr_cfa_pw <- as.matrix(loo_RLcr_cfa$pointwise[,3] )

fit_RLcr_2lr_cfa    <- readRDS("_outputs/fit_RLcumrew_2lr_cfa.RData")
lik_RLcr_2lr_cfa    <- extract_log_lik(fit_RLcr_2lr_cfa, parameter_name = "log_lik")
loo_RLcr_2lr_cfa    <- loo(lik_RLcr_2lr_cfa, cores = 3)
loo_RLcr_2lr_cfa_pw <- as.matrix(loo_RLcr_2lr_cfa$pointwise[,3] )

looMat = cbind(loo_RL_pw, loo_RLnc_pw, loo_RLnc_2lr_pw, loo_RLnc_cfa_pw, loo_RLnc_2lr_cfa_pw,
               loo_RLcoh_pw, loo_RLcoh_2lr_pw, loo_RLcoh_cfa_pw, loo_RLcoh_2lr_cfa_pw,
               loo_RLcr_pw, loo_RLcr_2lr_pw, loo_RLcr_cfa_pw, loo_RLcr_2lr_cfa_pw)
save(looMat, file = '_outputs/looicMat.RData')

#### plot ####
load(file = '_outputs/looicMat.RData')

looDF <- melt(looMat)
colnames(looDF) <- c('subID', 'modelName', 'LOOIC')
looDF$modelName <- c(rep('RL', 129), rep('RLnc',129), rep('RLnc_2lr',129),rep('RLnc_cfa',129),rep('RLnc_2lr_cfa',129),
                     rep('RLcoh',129),rep('RLcoh_2lr',129),rep('RLcoh_cfa',129),rep('RLcoh_2lr_cfa',129), 
                     rep('RLcumrew',129),rep('RLcumrew_2lr',129),rep('RLcumrew_cfa',129),rep('RLcumrew_2lr_cfa',129))
modelName <- c('RL','RLnc','RLnc_2lr','RLnc_cfa','RLnc_2lr_cfa',
               'RLcoh','RLcoh_2lr','RLcoh_cfa','RLcoh_2lr_cfa',
               'RLcumrew','RLcumrew_2lr','RLcumrew_cfa','RLcumrew_2lr_cfa')
looDF$modelName <- factor(looDF$modelName, levels = modelName)
theme_set(theme_gray(base_size = 20))

ggplot(data=looDF, aes(modelName, LOOIC)) + geom_boxplot(aes(fill = modelName)) 


