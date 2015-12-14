cal_prob_v1 <- function(dataList) {
  ## this function is part of the RLbeta_alt3 model
  ## it computes the action selection probability based on the learning rate and temperature 
  ## obtained from 'RevLearn_RLbeta_alt3_p1'
  ## the action probability is then plugged into 'RevLearn_RLbeata_alt3_p2'
  
  L <- list()
  ns <- dataList$nSubjects
  nt <- dataList$nTrials
  choice2 <- dataList$choice2
  wOthers <- dataList$wOthers
  otherChoice2 <- dataList$otherChoice2
  otherReward  <- dataList$otherReward
  
  load('_outputs/param_beta_alt3_v1.RData')
  lr <- param_beta_alt3[1:129,1]
  tau <- param_beta_alt3[130:258,1]
  
  pe   <- array(0,dim=c(nt,4))
  v    <- array(0,dim=c(nt+1,4,2))
  prob <- array(0,dim=c(ns,nt,4,2))
  prob_sC2  <- array(0,dim=c(ns,nt,4))
  prob_oC2  <- array(0,dim=c(ns,nt,4))
  wProb_sC2 <- array(0,dim=c(ns,nt,4))
  wProb_oC2 <- array(0,dim=c(ns,nt,4))
  
  for (s in 1:ns) {
    for (t in 1:nt) {
      for (o in 1:4) {
        prob[s,t,o,1] <- 1 / (1 + exp(tau[s] * (v[t,o,2]-v[t,o,1]) ))
        prob[s,t,o,2] <- 1 / (1 + exp(tau[s] * (v[t,o,1]-v[t,o,2]) ))
        pe[t,o]       <- otherReward[s,t,o] - v[t,o,otherChoice2[s,t,o]];
        v[t+1,o,] <- v[t,o,];
        v[t+1,o,otherChoice2[s,t,o]] <- v[t,o,otherChoice2[s,t,o]] + lr[s] * pe[t,o]; 
        
        ## calculated weighted action prob as the same/oppose to MY 2nd choice
        prob_sC2[s,t,o] <- prob[s,t,o,choice2[s,t]]
      }
      prob_oC2[s,t,]  <- 1 - prob_sC2[s,t,]
      wProb_sC2[s,t,] <- wOthers[s,t,] * prob_sC2[s,t,]
      wProb_oC2[s,t,] <- wOthers[s,t,] * prob_oC2[s,t,]
    } # trial loop
  }   # subject loop
  
  L$wProb_sC2 <- wProb_sC2
  L$wProb_oC2 <- wProb_oC2
  return(L)

} # fucntion
# end of function