run_model <- function(modelStr, test = TRUE, fitObj = NA, adapt = 0.8) {
  
  library(rstan); library(parallel); library(loo)
  L <- list()
    
  #### prepare data #### ===========================================================================
  dataList <- prep_data(modelStr)
    
  #### preparation for running stan #### ============================================================
  # model string in a separate .stan file
  modelFile <- paste0("_scripts/",modelStr,".stan")

  # setup up Stan configuration
  if (test == TRUE) {
    options(mc.cores = 1)
    nSamples <- 4
    nChains  <- 1 
    nBurnin  <- 0
    nThin    <- 1
  } else {
    options(mc.cores = 4)
    nSamples <- 4000#2000
    nChains  <- 4 
    nBurnin  <- floor(nSamples/2)
    nThin    <- 1#1
  }
  
  # parameter of interest (this could save both memory and space)
  poi <- create_pois(modelStr)

  #### run stan ####  ==============================================================================
  cat("Estimating", modelStr, "model... \n")
  startTime = Sys.time(); print(startTime)
  cat("Calling", nChains, "simulations in Stan... \n")
  rstan_options(auto_write = TRUE)
  
  stanfit <- stan(modelFile,
                fit     = fitObj,
                data    = dataList,
                pars    = poi,
                chains  = nChains,
                iter    = nSamples,
                warmup  = nBurnin,
                thin    = nThin,
                init    = "random",
                #seed    = 1581381385, # for RLbeta_alt3_p2_v2, LOOIV 8217
                #seed    = 295171225, # for RLbeta_alt2_v1, LOOIC 8216
                seed     = 1450154626, # for alt4_
                #seed     = 1136699103, # alt4_, seems even lower looic
                control = list(adapt_delta = adapt),
                verbose = FALSE)
  
  cat("Finishing", modelStr, "model simulation ... \n")
  endTime = Sys.time(); print(endTime)  
  cat("It took",as.character.Date(endTime - startTime), "\n")
  
  L$data <- dataList
  L$fit  <- stanfit
  
  return(L)
}  # function run_model()


#### nested functions #### ===========================================================================

prep_data <- function(modelstr){
  
  load("_data/sit_reversal_betnoshow_129.rdata")
  dataList <- list()
  sz <- dim(mydata)
  nt <- sz[1]; ns <- sz[3]
  dataList$nSubjects <- ns; dataList$nTrials <- nt
  
  choice1 <- array(0,dim = c(ns,nt)); choice2 <- array(0,dim = c(ns,nt)); reward <- array(0,dim = c(ns,nt))
  choice1 <- t(mydata[,3,])   # 1 OR  2
  choice2 <- t(mydata[,10,])  # 1 OR  2
  reward  <- t(mydata[,14,])  # 1 OR -1
  dataList$choice1 <- choice1
  dataList$choice2 <- choice2
  dataList$reward  <- reward
  
  if ( substr(modelstr,1,14) == "RevLearn_RLcoh" ) {
    
    chswtch <- array(0,dim = c(ns,nt))
    bet1    <- array(0,dim = c(ns,nt)); bet2    <- array(0,dim = c(ns,nt))
    with    <- array(0,dim = c(ns,nt)); against <- array(0,dim = c(ns,nt))
    my1     <- 0; other1  <- c(0,0,0,0)
    
    chswtch <- t(mydata[,5,])
    bet1    <- t(mydata[,13,]); bet2    <- t(mydata[,19,])
    
    for (s in 1:ns) {
      for (tr in 1:nt){
        my1 <- mydata[tr,3,s]; other1 <- mydata[tr,6:9,s]
        with[s,tr]    <- length(which(other1==my1)) /4  # count of with, either 1, 2, 3, or 4, divided by 4
        against[s,tr] <- length(which(other1!=my1)) /4  # count of against, either 1, 2, 3, or 4, divided by 4 
      }
    }
    
    dataList$chswtch <- chswtch
    dataList$bet1    <- bet1;    dataList$bet2    <- bet2
    dataList$with    <- with;    dataList$against <- against
    
  } else if ( substr(modelstr,1,20) == 'RevLearn_RLbeta_alt1' || substr(modelstr,1,20) == 'RevLearn_RLbeta_alt2' || 
              substr(modelstr,1,20) == 'RevLearn_RLbeta_alt3' || substr(modelstr,1,20) == 'RevLearn_RLbeta_alt4' ||
              substr(modelstr,1,20) == 'RevLearn_RLbeta_alt5') { 
    
    chswtch <- array(0,dim = c(ns,nt))
    bet1    <- array(0,dim = c(ns,nt)); bet2    <- array(0,dim = c(ns,nt))
    with    <- array(0,dim = c(ns,nt)); against <- array(0,dim = c(ns,nt))
    my1     <- 0; other1  <- c(0,0,0,0)
    otherChoice1 <- array(0,dim = c(ns,nt,4)); 
    otherChoice2 <- array(0,dim = c(ns,nt,4));
    otherReward  <- array(0,dim = c(ns,nt,4))
    pref         <- array(0,dim = c(ns,nt,4)); 
    wOthers      <- array(0,dim = c(ns,nt,4)) # others' weight [.75 .5 .25 .25]
    wOthers_one  <- array(0,dim = c(ns,nt,4)) # others' weight [ 3 2 1 1] / 7
    wghtValue    <- array(0,dim = c(ns,nt,2)) # others' value based on weight
    cfsC2        <- array(0,dim = c(ns,nt,4)) # cumulative-window frequency, same as my C2 
    cfoC2        <- array(0,dim = c(ns,nt,4)) # cumulative-window frequency, opposite to my C2
    wgtWith      <- array(0,dim = c(ns,nt)); wgtWith_one <- array(0,dim = c(ns,nt))
    wgtAgst      <- array(0,dim = c(ns,nt)); wgtAgst_one <- array(0,dim = c(ns,nt))
    otherCumAcc  <- array(0,dim = c(ns,nt,4));  # cumulative accuracy according to reward probability
    
    chswtch <- t(mydata[,5,])
    bet1    <- t(mydata[,13,]); bet2    <- t(mydata[,19,])
    wgtWith <- t(mydata[,93,]); wgtAgst <- t(mydata[,94,])
    wgtWith_one <- t(mydata[,99,])
    wgtAgst_one <- t(mydata[,100,])
    
    for (s in 1:ns) {
      otherChoice1[s,,] <- mydata[,6:9,s]
      otherChoice2[s,,] <- mydata[,55:58,s]
      otherReward[s,,]  <- mydata[,24:27,s]
      pref[s,,]         <- mydata[,47:50,s]
      wOthers[s,,]      <- mydata[,51:54,s]
      wOthers_one[s,,]  <- mydata[,95:98,s]
      wghtValue[s,,]    <- mydata[,59:60,s]
      cfsC2[s,,]        <- mydata[,61:64,s]
      cfoC2[s,,]        <- mydata[,65:68,s]
      otherCumAcc[s,,]  <- mydata[,129:132,s]
      
      for (t in 1:nt){
        my1 <- mydata[t,3,s]; other1 <- mydata[t,6:9,s]
        with[s,t]    <- length(which(other1==my1)) /4
        against[s,t] <- length(which(other1!=my1)) /4
      }
    }
    dataList$chswtch      <- chswtch;
    dataList$otherChoice1 <- otherChoice1
    dataList$otherChoice2 <- otherChoice2
    dataList$otherReward  <- otherReward
    dataList$wghtValue    <- wghtValue
    dataList$bet1    <- bet1;    dataList$bet2    <- bet2
    dataList$with    <- with;    dataList$against <- against
    dataList$pref    <- pref;    dataList$wOthers <- wOthers
    dataList$cfsC2   <- cfsC2;   dataList$cfoC2   <- cfoC2;
    dataList$wgtWith <- wgtWith; dataList$wgtAgst <- wgtAgst
    dataList$wOthers_one <- wOthers_one
    dataList$wgtWith_one <- wgtWith_one
    dataList$wgtAgst_one <- wgtAgst_one
    dataList$otherCumAcc <- otherCumAcc
    
    if ( substr(modelstr,1,20) == 'RevLearn_RLbeta_alt2' ) {
      wProb_sC2 <- array(0,dim = c(ns,nt,4)) 
      wProb_oC2 <- array(0,dim = c(ns,nt,4))
      wProb_sC2_med <- array(0,dim = c(ns,nt,4)) 
      wProb_oC2_med <- array(0,dim = c(ns,nt,4))
      
      for (s in 1:ns) {
        wProb_sC2[s,,] <- mydata[,81:84,s]
        wProb_oC2[s,,] <- mydata[,85:88,s]
        wProb_sC2_med[s,,] <- mydata[,109:112,s]
        wProb_oC2_med[s,,] <- mydata[,113:116,s]
      }
      dataList$wProb_sC2 <- wProb_sC2
      dataList$wProb_oC2 <- wProb_oC2
      dataList$wProb_sC2_med <- wProb_sC2_med
      dataList$wProb_oC2_med <- wProb_oC2_med
      
    } else if (modelstr == "RevLearn_RLbeta_alt3_p2_v1" || modelstr == "RevLearn_RLbeta_alt3_p2_v1_w") {
      L <- cal_prob_v1(dataList)
      dataList$wProb_sC2 <- L$wProb_sC2
      dataList$wProb_oC2 <- L$wProb_oC2
    } else if (modelstr == "RevLearn_RLbeta_alt3_p2_v2") {
      L <- cal_prob_v2(dataList)
      dataList$wProb_sC2 <- L$wProb_sC2
      dataList$wProb_oC2 <- L$wProb_oC2
    } else if ( substr(modelstr,1,20) == 'RevLearn_RLbeta_alt4' || substr(modelstr,1,20) == 'RevLearn_RLbeta_alt5' ) {
      otherReward2 <- array(0,dim = c(ns,nt,4))
      otherReward2_actural <- array(0,dim = c(ns,nt,4))
      otherWith2   <- array(0,dim = c(ns,nt,4))
      for (s in 1:ns) {
        otherReward2[s,,]  <- mydata[,24:27,s]
        otherWith2[s,,]    <- mydata[,89:92,s]  # otherChoice2 == myChoice2, with(1) or against(0)
      }
      otherReward2_actural <- otherReward2
      
      if (modelstr != "RevLearn_RLbeta_alt4_c_w_v12_1lr" && modelstr != "RevLearn_RLbeta_alt4_c_w_v21_1lr" &&
          modelstr != "RevLearn_RLbeta_alt4_c_w_v25_1lr" && modelstr != "RevLearn_RLbeta_alt4_c_w_v27_1lr") {
          otherReward2[otherReward2 == -1] = 0
      }
      dataList$otherReward2  <- otherReward2  # [0 1] OR [-1 1]
      dataList$otherReward2_actural  <- otherReward2_actural   ## always [-1 1]
      dataList$otherWith2    <- otherWith2
    } 
    
    
  } else if ( substr(modelstr,1,17) == "RevLearn_RLcumrew" ) {
    
    otherChoice1 <- array(0,dim = c(ns,nt,4))
    otherReward  <- array(0,dim = c(ns,nt,4))
    otherWith    <- array(0,dim = c(ns,nt,4))
    
    for (s in 1:ns) {
      otherChoice1[s,,] <- mydata[,6:9,s]
      otherReward[s,,]  <- mydata[,24:27,s]
      otherWith[s,,]    <- mydata[,69:72,s]  # otherChoice1 == myChoice1, with(1) or against(0)
    }
    otherReward[otherReward == -1] = 0
    dataList$otherChoice1 <- otherChoice1
    dataList$otherReward  <- otherReward
    dataList$otherWith    <- otherWith
  }
  
  return(dataList)
} # function


# ------------------------------------------------------------------------------------------------------
create_pois <- function(model){
  pois <- list()
  
  if (model == "RevLearn_RL" || model == "RevLearn_RLnc"){
    pois <- c("lr_mu", "tau_mu", 
              "lr_sd", "tau_sd",
              "lr", "tau", 
              "c_rep",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLnc_2lr") {
    pois <- c("lr1_mu", "lr2_mu", "tau_mu",  
              "lr1_sd", "lr2_sd", "tau_sd", 
              "lr1", "lr2", "tau",
              "c_rep",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLnc_cfa") {
    pois <- c("lr_mu", "tau_mu",  "cfa_mu", 
              "lr_sd", "tau_sd", "cfa_sd",
              "lr", "tau", "cfa",
              "c_rep",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLnc_2lr_cfa") {
    pois <- c("lr1_mu", "lr2_mu", "tau_mu", "cfa_mu",
              "lr1_sd", "lr2_sd", "tau_sd", "cfa_sd",
              "lr1", "lr2", "tau", "cfa",
              "c_rep",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLcoh"){
    pois <- c("lr_mu", "tau_mu", "coha_mu", "cohw_mu", 
              "lr_sd", "tau_sd", "coha_sd", "cohw_sd",
              "lr", "tau", "coha", "cohw",
              "c_rep",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLcoh_2lr"){
    pois <- c("lr1_mu", "lr2_mu", "tau_mu", "coha_mu", "cohw_mu", 
              "lr1_sd", "lr2_sd", "tau_sd", "coha_sd", "cohw_sd",
              "lr1", "lr2", "tau", "coha", "cohw",
              "c_rep",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLcoh_cfa") {
    pois <- c("lr_mu", "tau_mu", "coha_mu", "cohw_mu", "cfa_mu",
              "lr_sd", "tau_sd", "coha_sd", "cohw_sd", "cfa_sd",
              "lr", "tau", "coha", "cohw", "cfa",
              "c_rep",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLcoh_2lr_cfa") {
    pois <- c("lr1_mu", "lr2_mu", "tau_mu", "coha_mu", "cohw_mu", "cfa_mu",
              "lr1_sd", "lr2_sd", "tau_sd", "coha_sd", "cohw_sd", "cfa_sd",
              "lr1", "lr2", "tau", "coha", "cohw", "cfa", 
              "c_rep",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLcoh_modvalue" || model == "RevLearn_RLcoh_modprob_tempin" || model == "RevLearn_RLcoh_modprob_tempout"){
    pois <- c("lr_mu", "tau_mu", "coha_mu", "cohw_mu", 
              "lr_sd", "tau_sd", "coha_sd", "cohw_sd",
              "lr", "tau", "coha", "cohw",
              "log_lik1", 
              "log_lik2", "lp__")
  } else if (model == "RevLearn_RLcoh_2lr2t_modvalue" || model == "RevLearn_RLcoh_2lr2t_modprob_tempin" || 
             model == "RevLearn_RLcoh_2lr2t_modprob_tempout" || model == "RevLearn_RLcoh_2lr2t_modprob_tempin_bern") {   
    pois <- c("lr1_mu", "tau1_mu", "lr2_mu", "tau2_mu", "coha_mu", "cohw_mu", 
              "lr1_sd", "tau1_sd", "lr2_sd", "tau2_sd", "coha_sd", "cohw_sd",
              "lr1", "tau1", "lr2", "tau2", "coha", "cohw",
              "log_lik1", 
              "log_lik2", "lp__")
  } else if (model == "RevLearn_RLcumrew" ) {
    pois <- c("lr_mu", "tau_mu", "disc_mu", "cra_mu", "crw_mu", 
              "lr_sd", "tau_sd", "disc_sd", "cra_sd", "crw_sd",
              "lr", "tau", "disc", "cra", "crw",
              "c_rep",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLcumrew_2lr" ) {
    pois <- c("lr1_mu", "lr2_mu", "tau_mu", "disc_mu", "cra_mu", "crw_mu", 
              "lr1_sd", "lr2_sd", "tau_sd", "disc_sd", "cra_sd", "crw_sd",
              "lr1", "lr2", "tau", "disc", "cra", "crw",
              "c_rep",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLcumrew_cfa" ) {
    pois <- c("lr_mu", "tau_mu", "disc_mu", "cra_mu", "crw_mu", "cfa_mu",
              "lr_sd", "tau_sd", "disc_sd", "cra_sd", "crw_sd", "cfa_sd",
              "lr", "tau", "disc", "cra", "crw", "cfa",
              "c_rep",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLcumrew_2lr_cfa" ) {
    pois <- c("lr1_mu", "lr2_mu", "tau_mu", "disc_mu", "cra_mu", "crw_mu", "cfa_mu",
              "lr1_sd", "lr2_sd", "tau_sd", "disc_sd", "cra_sd", "crw_sd", "cfa_sd",
              "lr1", "lr2", "tau", "disc", "cra", "crw", "cfa",
              "c_rep",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLbeta_alt1_bc") {
    pois <- c("lr_mu", "thrs_mu", "beta_mu",
              "lr_sd", "thrs_sd", "beta_sd",
              "lr", "thrs", "beta",
              "log_likc1", "log_likc2", "log_likb1", "log_likb2", "lp__")
  } else if ( substr(model,1,22) == 'RevLearn_RLbeta_alt1_c' || substr(model,1,25) == 'RevLearn_RLbeta_alt2_c_v2' || 
              substr(model,1,23) == 'RevLearn_RLbeta_alt3_p2' ) {
    pois <- c("lr_mu", "beta_mu",
              "lr_sd", "beta_sd",
              "lr", "beta",
              "c_rep",
              "log_likc1", "log_likc2", "lp__")
  } else if (model == "RevLearn_RLbeta_alt2_bc") {
    pois <- c("lr_mu", "thrs_mu", "evid_wght_mu", "beta_mu",
              "lr_sd", "thrs_sd", "evid_wght_sd", "beta_sd",
              "lr", "thrs", "evid_wght", "beta",
              "log_likc1", "log_likc2", "log_likb1", "log_likb2", "lp__")
  } else if ( substr(model,1,25) == 'RevLearn_RLbeta_alt2_c_v1' ) {
    pois <- c("lr_mu", "evidW_mu", "beta_mu",
              "lr_sd", "evidW_sd", "beta_sd",
              "lr", "evidW", "beta",
              "c_rep",
              "log_likc1", "log_likc2", "lp__")
  } else if (model == "RevLearn_RLbeta_alt3_p1_v1" || model == "RevLearn_RLbeta_alt3_p1_v2") {
    pois <- c("lr_mu", "tau_mu", 
              "lr_sd", "tau_sd",
              "lr", "tau", 
              "lp__")
  } else if ( model == 'RevLearn_RLbeta_alt4_c_w_v10_1lr' ) {
    pois <- c("lr_mu", "beta_mu", "disc_mu", "cfa_mu",
              "lr_sd", "beta_sd", "disc_sd", "cfa_sd",
              "lr", "beta", "disc", "cfa",
              "c_rep",
              "log_likc1", "log_likc2", "lp__")
    
  } else if ( model == 'RevLearn_RLbeta_alt4_c_w_v15_1lr' ) {
    pois <- c("lr", "beta", "disc",
              "c_rep",
              "log_likc1", "log_likc2", "lp__")
    
  } else if ( model == 'RevLearn_RLbeta_alt4_c_w_v23_1lr' ) {
    pois <- c("lr_mu", "beta_mu", "disc_mu", "tau_mu",
              "lr_sd", "beta_sd", "disc_sd", "tau_sd",
              "lr", "beta", "disc", "tau",
              "c_rep",
              "log_likc1", "log_likc2", "lp__")
  } else if ( substr(model,1,20) == 'RevLearn_RLbeta_alt4' || substr(model,1,20) == 'RevLearn_RLbeta_alt5' ) {
    pois <- c("lr_mu", "beta_mu", "disc_mu",
              "lr_sd", "beta_sd", "disc_sd",
              "lr", "beta", "disc",
              "c_rep",
              "log_likc1", "log_likc2", "lp__")
  }
  
  return(pois)
} # function

#### end of function ####