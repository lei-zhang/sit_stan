fit_RL    <- readRDS("_outputs/fit_RL.RData")
loo_RL    <- extract_looic(fit_RL,NULL)$looic
aic_RL    <- extract_otheric(fit_RL,NULL, np = 2)$myAIC
bic_RL    <- extract_otheric(fit_RL,NULL, np = 2)$myBIC
waic_RL   <- extract_otheric(fit_RL,NULL, np = 2)$myWAIC$waic
# rm(fit_RL)

fit_RLnc  <- readRDS("_outputs/fit_RLnc.RData")
loo_RLnc  <- extract_looic(fit_RLnc,NULL)$looic
aic_RLnc    <- extract_otheric(fit_RLnc,NULL, np = 2)$myAIC
bic_RLnc    <- extract_otheric(fit_RLnc,NULL, np = 2)$myBIC
waic_RLnc   <- extract_otheric(fit_RLnc,NULL, np = 2)$myWAIC$waic
#rm(fit_RLnc)

fit_RLcoh <- readRDS("_outputs/fit_RLcoh.RData")
loo_RLcoh <- extract_looic(fit_RLcoh,NULL)$looic
aic_RLcoh    <- extract_otheric(fit_RLcoh,NULL, np = 4)$myAIC
bic_RLcoh    <- extract_otheric(fit_RLcoh,NULL, np = 4)$myBIC
waic_RLcoh   <- extract_otheric(fit_RLcoh,NULL, np = 4)$myWAIC$waic
# rm(fit_RLcoh)

# fit_RLcr  <- readRDS("_outputs/fit_RLcumrew.RData")
# loo_RLcr <- extract_looic(fit_RLcr,NULL)$looic
# rm(fit_RLcr)

fit_RLbeta_a1   <- readRDS("_outputs/fit_RLbeta_alt1_c_w_080.RData")
loo_RLbeta_a1   <- extract_looic(fit_RLbeta_a1,core=1)$looic
aic_RLbeta_a1  <- extract_otheric(fit_RLbeta_a1,np = 7)$myAIC
bic_RLbeta_a1  <- extract_otheric(fit_RLbeta_a1,np = 7)$myBIC
waic_RLbeta_a1 <- extract_otheric(fit_RLbeta_a1,np = 7)$myWAIC$waic
#rm(fit_RLbeta_a1)

fit_RLbeta_a2   <- readRDS("_outputs/fit_RLbeta_alt2_c_v2_w_1lr.RData")
loo_RLbeta_a2   <- extract_looic(fit_RLbeta_a2,core=1)$looic
aic_RLbeta_a2  <- extract_otheric(fit_RLbeta_a2,np = 7)$myAIC
bic_RLbeta_a2  <- extract_otheric(fit_RLbeta_a2,np = 7)$myBIC
waic_RLbeta_a2 <- extract_otheric(fit_RLbeta_a2,np = 7)$myWAIC$waic
#rm(fit_RLbeta_a2)

fit_RLbeta_a3   <- readRDS("_outputs/fit_RLbeta_alt3_p2_v1_w_080.RData")
loo_RLbeta_a3   <- extract_looic(fit_RLbeta_a3,core=1)$looic
aic_RLbeta_a3  <- extract_otheric(fit_RLbeta_a3,np = 9)$myAIC
bic_RLbeta_a3  <- extract_otheric(fit_RLbeta_a3,np = 9)$myBIC
waic_RLbeta_a3 <- extract_otheric(fit_RLbeta_a3,np = 9)$myWAIC$waic
#rm(fit_RLbeta_a3)

fit_RLbeta_a4   <- readRDS("_outputs/fit_RLbeta_alt4_c_w_v6.RData")
loo_RLbeta_a4   <- extract_looic(fit_RLbeta_a4,core=1)$looic
aic_RLbeta_a4  <- extract_otheric(fit_RLbeta_a4,np = 8)$myAIC
bic_RLbeta_a4  <- extract_otheric(fit_RLbeta_a4,np = 8)$myBIC
waic_RLbeta_a4 <- extract_otheric(fit_RLbeta_a4,np = 8)$myWAIC$waic
# rm(fit_RLbeta_a4)

AIC_vec   <- c(aic_RL, aic_RLnc, aic_RLcoh, aic_RLbeta_a1, aic_RLbeta_a2, aic_RLbeta_a3, aic_RLbeta_a4)
BIC_vec   <- c(bic_RL, bic_RLnc, bic_RLcoh, bic_RLbeta_a1, bic_RLbeta_a2, bic_RLbeta_a3, bic_RLbeta_a4)
WAIC_vec  <- c(waic_RL, waic_RLnc, waic_RLcoh, waic_RLbeta_a1, waic_RLbeta_a2, waic_RLbeta_a3, waic_RLbeta_a4)
LOOIC_vec <- c(loo_RL, loo_RLnc, loo_RLcoh, loo_RLbeta_a1, loo_RLbeta_a2, loo_RLbeta_a3, loo_RLbeta_a4)

IC_mat <- cbind(AIC_vec, BIC_vec, WAIC_vec, LOOIC_vec)

save(IC_mat, file = '_outputs/ICMat.RData')