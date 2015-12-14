data {
  int<lower=1> nSubjects;                                    // number of subjects
  int<lower=1> nTrials;                                      // number of trials 
  int<lower=1,upper=2> choice1[nSubjects,nTrials];           // 1st choices, 1 or 2
  int<lower=1,upper=2> choice2[nSubjects,nTrials];           // 2nd choices, 1 or 2
  int<lower=0,upper=1> chswtch[nSubjects,nTrials];           // choice switch, 0 or 1
  real<lower=-1,upper=1> reward[nSubjects,nTrials];          // outcome, 1 or -1
  real<lower=0,upper=1> with[nSubjects,nTrials];             // No. of with
  real<lower=0,upper=1> against[nSubjects,nTrials];          // No. of against
  real<lower=0.25,upper=1> wOthers[nSubjects,nTrials,4];     // others' weight based on pref
  real<lower=0,upper=5> cfsC2[nSubjects,nTrials,4];          // cum-window freq, same as my C2
  real<lower=0,upper=5> cfoC2[nSubjects,nTrials,4];          // cum-window freq, opposite to my C2
}

transformed data {
  vector[2] initV;  // initial values for V
  int B;            // number of predictors
  
  initV <- rep_vector(0.0,2); 
  B     <- 6;
}

parameters {
  // group-level parameters
  real lr_mu_pr;    
  real evidW_mu_pr;
  vector[B] beta_mu;
  
  real<lower=0> lr_sd;
  real<lower=0> evidW_sd;
  vector<lower=0>[B] beta_sd;
  
  // subject-level raw parameters, follows norm(0,1), for later Matt Trick
  vector[nSubjects] lr_raw;
  vector[nSubjects] evidW_raw;
  vector[nSubjects] beta_raw[B];      // dim: [B nSubjects]
}

transformed parameters {
  // subject-level parameters
  vector<lower=0,upper=1>[nSubjects] lr;
  vector<lower=0,upper=1>[nSubjects] evidW;
  vector[nSubjects] beta[B];
  
  // Matt Trick, note that the input of Phi_approx must be 'real' rather than 'vector'
  for (s in 1:nSubjects) {
    lr[s]    <- Phi_approx( lr_mu_pr + lr_sd * lr_raw[s] );
    evidW[s] <- Phi_approx( evidW_mu_pr + evidW_sd * evidW_raw[s] );
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
  vector[4] prob_sC2;  // others' prob of choosing the same as my C2
  vector[4] prob_oC2;  // others' prob of choosing the opposite as my C2
  vector[4] wProb_sC2; // weighted prob based on preference' weight, same as my C2
  vector[4] wProb_oC2; // weighted prob based on preference' weight, opposite to my C2
  
  // hyperparameters
  lr_mu_pr    ~ normal(0,1);
  evidW_mu_pr ~ normal(0,1);
  beta_mu     ~ normal(0,1);
  
  lr_sd    ~ cauchy(0,5);
  evidW_sd ~ cauchy(0,5);
  beta_sd  ~ cauchy(0,5);
  
  // Matt Trick
  lr_raw    ~ normal(0,1);
  evidW_raw ~ normal(0,1);
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
      chswtch[s,t] ~ bernoulli_logit( valfun2 );
      
      // my prediction error
      pe[t]   <-  reward[s,t] - myValue[t,choice2[s,t]];
      penc[t] <- -reward[s,t] - myValue[t,3-choice2[s,t]];
      
      // update my value
      myValue[t+1,choice2[s,t]]   <- myValue[t,choice2[s,t]]   + lr[s] * pe[t];
      myValue[t+1,3-choice2[s,t]] <- myValue[t,3-choice2[s,t]] + lr[s] * penc[t];
      
      // use choice preference from others to update the others' value
      for (o in 1:4) {
        // [note] the [a b] parameter of beta_cdf must be positive, and count from beta(1,1)
        prob_sC2[o]  <- beta_cdf(0.5, evidW[s]*cfoC2[s,t,o] + 1, evidW[s]*cfsC2[s,t,o] + 1 );        
      }
      prob_oC2  <- 1 - prob_sC2;
      wProb_sC2 <- to_vector(wOthers[s,t]) .* prob_sC2;
      wProb_oC2 <- to_vector(wOthers[s,t]) .* prob_oC2;

      otherValue[t+1,choice2[s,t]]   <- sum( wProb_sC2 );
      otherValue[t+1,3-choice2[s,t]] <- sum( wProb_oC2 );
    }  // trial loop
  }    // subject loop
}

generated quantities {
  real<lower=0,upper=1> lr_mu;
  real<lower=0,upper=1> evidW_mu;
  
  real log_likc1[nSubjects];
  real log_likc2[nSubjects];
  vector[2] myValue2[nTrials+1];
  vector[2] otherValue2[nTrials+1];
  vector[nTrials] pe2;
  vector[nTrials] penc2;
  real valdiff_gen;
  vector[2] valfun1_gen;
  real valfun2_gen;
  vector[4] prob_sC2_gen;
  vector[4] prob_oC2_gen;
  vector[4] wProb_sC2_gen;
  vector[4] wProb_oC2_gen;
  int<lower=0,upper=1> c_rep[nSubjects, nTrials];
  
  lr_mu    <- Phi_approx(lr_mu_pr);
  evidW_mu <- Phi_approx(evidW_mu_pr);
  
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
      myValue2[t+1,choice2[s,t]]   <- myValue2[t,choice2[s,t]]   + lr[s] * pe2[t];
      myValue2[t+1,3-choice2[s,t]] <- myValue2[t,3-choice2[s,t]] + lr[s] * penc2[t];

      for (o in 1:4) {
        prob_sC2_gen[o]  <- beta_cdf(0.5, evidW[s]*cfoC2[s,t,o] + 1, evidW[s]*cfsC2[s,t,o] + 1 );
      }
      prob_oC2_gen  <- 1 - prob_sC2_gen;
      wProb_sC2_gen <- to_vector(wOthers[s,t]) .* prob_sC2_gen;
      wProb_oC2_gen <- to_vector(wOthers[s,t]) .* prob_oC2_gen;

      otherValue2[t+1,choice2[s,t]]   <- sum( wProb_sC2_gen );
      otherValue2[t+1,3-choice2[s,t]] <- sum( wProb_oC2_gen );
    }  // trial loop
  }    // subject loop
}