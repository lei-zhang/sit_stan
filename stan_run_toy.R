run_model_toy <- function(modelStr, test = TRUE, fitObj = NA, adapt = 0.8) {
  
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
    nSamples <- 2000#2000
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
  mydata   <- mydata[,,62:63]
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
  
  if ( substr(modelstr,1,20) == 'RevLearn_RLbeta_alt4' || substr(modelstr,1,20) == 'RevLearn_RLbeta_alt5' ) {
    otherReward2 <- array(0,dim = c(ns,nt,4))
    otherWith2   <- array(0,dim = c(ns,nt,4))
    for (s in 1:ns) {
      otherReward2[s,,]  <- mydata[,24:27,s]
      otherWith2[s,,]    <- mydata[,89:92,s]  # otherChoice2 == myChoice2, with(1) or against(0)
    }
    if (modelstr != "RevLearn_RLbeta_alt4_c_w_v12_1lr" && modelstr != "RevLearn_RLbeta_alt4_c_w_v21_1lr" &&
        modelstr != "RevLearn_RLbeta_alt4_c_w_v25_1lr" && modelstr != "RevLearn_RLbeta_alt4_c_w_v27_1lr") {
      otherReward2[otherReward2 == -1] = 0
    }
    dataList$otherReward2  <- otherReward2
    dataList$otherWith2    <- otherWith2
  } 
  
  return(dataList)
} # function


# ------------------------------------------------------------------------------------------------------
create_pois <- function(model){
  pois <- list()
  
  if ( model == 'RevLearn_RLbeta_alt4_c_w_v10_1lr' ) {
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