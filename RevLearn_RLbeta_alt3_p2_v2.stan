data {
  int<lower=1> nSubjects;                                // number of subjects
  int<lower=1> nTrials;                                  // number of trials 
  int<lower=1,upper=2> choice1[nSubjects,nTrials];       // 1st choices, 1 or 2
  int<lower=1,upper=2> choice2[nSubjects,nTrials];       // 2nd choices, 1 or 2
  int<lower=0,upper=1> chswtch[nSubjects,nTrials];       // choice switch, 0 or 1
  real<lower=-1,upper=1> reward[nSubjects,nTrials];      // outcome, 1 or -1
  real<lower=0,upper=4> with[nSubjects,nTrials];         // No. of with
  real<lower=0,upper=4> against[nSubjects,nTrials];      // No. of against
  real<lower=0,upper=1> wProb_sC2[nSubjects,nTrials,4];  // weighted prob based on preference' weight, same as my C2
  real<lower=0,upper=1> wProb_oC2[nSubjects,nTrials,4];  // weighted prob based on preference' weight, opposite to my C2
}

transformed data {
  vector[2] initV;  // initial values for V
  int B;            // number of predictors

  initV <- rep_vector(0.0,2); 
  B     <- 6;         // number of predictors
}

parameters {
  // group-level parameters
  //vector[2] lr_mu_pr;
  real lr_mu_pr;
  vector[B] beta_mu;
  
  //vector<lower=0>[2] lr_sd;
  real<lower=0> lr_sd;
  vector<lower=0>[B] beta_sd;
  
  // subject-level raw parameters, follows norm(0,1), for later Matt Trick
  //vector[nSubjects] lr_raw[2];        // dim: [2 nSubjects]
  vector[nSubjects] lr_raw;
  vector[nSubjects] beta_raw[B];      // dim: [B nSubjects]
}

transformed parameters {
  // subject-level parameters
  //vector<lower=0,upper=1>[nSubjects] lr[2];
  vector<lower=0,upper=1>[nSubjects] lr;
  vector[nSubjects] beta[B];
  
  // Matt Trick, note that the input of Phi_approx must be 'real' rather than 'vector'
//   for (i in 1:2)           {
//     for (s in 1:nSubjects) {
//       lr[i,s]    <- Phi_approx( lr_mu_pr[i] + lr_sd[i] * lr_raw[i,s] );
//     }
//   }
  for (s in 1:nSubjects) {
    lr[s]  <- Phi_approx( lr_mu_pr + lr_sd * lr_raw[s] );
  }
  
  for (i in 1:B) {  // partial vectorization
    beta[i] <- beta_mu[i] + beta_sd[i] * beta_raw[i];
  }
}

model {
  // define the value and pe vectors
  vector[2] myValue[nTrials+1];
  vector[2] otherValue[nTrials+1];
  vector[nTrials] pe;
  vector[nTrials] penc;
  real valdiff;
  vector[2] valfun1;
  real valfun2;
  
  // hyperparameters
  lr_mu_pr    ~ normal(0,1);
  beta_mu     ~ normal(0,1);
  
  lr_sd    ~ cauchy(0,5);
  beta_sd  ~ cauchy(0,5);
  
  // Matt Trick
//   for (i in 1:2) {
//     lr_raw[i]    ~ normal(0,1);
//   }
  lr_raw   ~ normal(0,1);
  
  for (i in 1:B) {
    beta_raw[i]  ~ normal(0,1);
  }
  
  // subject loop and trial loop
  for (s in 1:nSubjects) {
    myValue[1]    <- initV;
    otherValue[1] <- initV;
    
    for (t in 1:nTrials) {
      valfun1 <- beta[1,s]*myValue[t] + beta[2,s]*otherValue[t];
      choice1[s,t] ~ categorical_logit( valfun1 );
      
      valdiff <- myValue[t,choice1[s,t]] - myValue[t,3-choice1[s,t]];
      valfun2 <- beta[3,s] + beta[4,s]*valdiff + beta[5,s]*with[s,t] + beta[6,s]*against[s,t];
      chswtch[s,t] ~ bernoulli_logit(valfun2);
      
      // my prediction error
      pe[t]   <-  reward[s,t] - myValue[t,choice2[s,t]];
      penc[t] <- -reward[s,t] - myValue[t,3-choice2[s,t]];
      
      // update my value
      //myValue[t+1,choice2[s,t]]   <- myValue[t,choice2[s,t]]   + lr[1,s] * pe[t];
      //myValue[t+1,3-choice2[s,t]] <- myValue[t,3-choice2[s,t]] + lr[2,s] * penc[t];
      myValue[t+1,choice2[s,t]]   <- myValue[t,choice2[s,t]]   + lr[s] * pe[t];
      myValue[t+1,3-choice2[s,t]] <- myValue[t,3-choice2[s,t]] + lr[s] * penc[t];
      
      // treat other co-players as a separated reinforcer to update the others' value
      otherValue[t+1,choice2[s,t]]   <- sum( wProb_sC2[s,t] );
      otherValue[t+1,3-choice2[s,t]] <- sum( wProb_oC2[s,t] );
    }  // trial loop
  }    // subject loop
}

generated quantities {
  //real<lower=0,upper=1> lr_mu[2];
  real<lower=0,upper=1> lr_mu; 
  
  real log_likc1[nSubjects];
  real log_likc2[nSubjects];
  vector[2] myValue2[nTrials+1];
  vector[2] otherValue2[nTrials+1];
  vector[nTrials] pe2;
  vector[nTrials] penc2;
  real valdiff_gen;
  vector[2] valfun1_gen;
  real valfun2_gen;
  int<lower=0,upper=1> c_rep[nSubjects, nTrials];
  
  //for (i in 1:2) {
  //  lr_mu[i]    <- Phi_approx(lr_mu_pr[i]);    
  //}
  lr_mu  <- Phi_approx(lr_mu_pr);
  
  for (s in 1:nSubjects) {
    myValue2[1]    <- initV;
    otherValue2[1] <- initV;
    log_likc1[s]   <- 0;
    log_likc2[s]   <- 0;
    
    for (t in 1:nTrials) {
      valfun1_gen  <- beta[1,s]*myValue2[t] + beta[2,s]*otherValue2[t];
      log_likc1[s] <- log_likc1[s] + categorical_logit_log(choice1[s,t], valfun1_gen);
      
      valdiff_gen  <- myValue2[t,choice1[s,t]] - myValue2[t,3-choice1[s,t]];
      valfun2_gen  <- beta[3,s] + beta[4,s]*valdiff_gen + beta[5,s]*with[s,t] + beta[6,s]*against[s,t];
      log_likc2[s] <- log_likc2[s] + bernoulli_logit_log(chswtch[s,t], valfun2_gen);

      c_rep[s,t]   <- bernoulli_rng( inv_logit(valfun2_gen) );
      
      pe2[t]   <-  reward[s,t] - myValue2[t,choice2[s,t]];
      penc2[t] <- -reward[s,t] - myValue2[t,3-choice2[s,t]];
      
      //myValue2[t+1,choice2[s,t]]   <- myValue2[t,choice2[s,t]]   + lr[1,s] * pe2[t];
      //myValue2[t+1,3-choice2[s,t]] <- myValue2[t,3-choice2[s,t]] + lr[2,s] * penc2[t];
      myValue2[t+1,choice2[s,t]]   <- myValue2[t,choice2[s,t]]   + lr[s] * pe2[t];
      myValue2[t+1,3-choice2[s,t]] <- myValue2[t,3-choice2[s,t]] + lr[s] * penc2[t];

      otherValue2[t+1,choice2[s,t]]   <- sum( wProb_sC2[s,t] );
      otherValue2[t+1,3-choice2[s,t]] <- sum( wProb_oC2[s,t] );
    }  // trial loop
  }    // subject loop
}