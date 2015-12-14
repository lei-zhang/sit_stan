fit_RL    <- readRDS("_outputs/fit_RL.RData")
loo_RL    <- extract_looic(fit_RL,NULL)$looic
#rm(fit_RL)

fit_RLnc  <- readRDS("_outputs/fit_RLnc.RData")
loo_RLnc  <- extract_looic(fit_RLnc,NULL)$looic
#rm(fit_RLnc)

fit_RLcoh <- readRDS("_outputs/fit_RLcoh.RData")
loo_RLcoh <- extract_looic(fit_RLcoh,NULL)$looic
#rm(fit_RLcoh)

fit_RLcr  <- readRDS("_outputs/fit_RLcumrew.RData")
loo_RLcr <- extract_looic(fit_RLcr,NULL)$looic
# rm(fit_RLcr)

fit_RLbeta_a1 <- readRDS("_outputs/fit_RLbeta_alt1_c_w_080.RData")
loo_RLbeta_a1 <- extract_looic(fit_RLbeta_a1,core=1)$looic
rm(fit_RLbeta_a1)

fit_RLbeta_a2 <- readRDS("_outputs/fit_RLbeta_alt2_c_v2_w_1lr.RData")
loo_RLbeta_a2 <- extract_looic(fit_RLbeta_a2,core=1)$looic
rm(fit_RLbeta_a2)

fit_RLbeta_a3 <- readRDS("_outputs/fit_RLbeta_alt3_p2_v1_w_080.RData")
loo_RLbeta_a3 <- extract_looic(fit_RLbeta_a3,core=1)$looic
rm(fit_RLbeta_a3)

fit_RLbeta_a4 <- readRDS("_outputs/fit_RLbeta_alt4_c_w_v6.RData")
loo_RLbeta_a4 <- extract_looic(fit_RLbeta_a4,core=1)$looic
rm(fit_RLbeta_a4)



# looMat = cbind(loo_RL_pw, loo_RLnc_pw, loo_RLnc_2lr_pw, loo_RLnc_cfa_pw, loo_RLnc_2lr_cfa_pw,
#                loo_RLcoh_pw, loo_RLcoh_2lr_pw, loo_RLcoh_cfa_pw, loo_RLcoh_2lr_cfa_pw,
#                loo_RLcr_pw, loo_RLcr_2lr_pw, loo_RLcr_cfa_pw, loo_RLcr_2lr_cfa_pw)
# save(looMat, file = '_outputs/looicMat.RData')