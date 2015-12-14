data {
  int<lower=1> nSubjects;                              // number of subjects
  int<lower=1> nTrials;                                // number of trials 
  int<lower=1,upper=2> choice1[nSubjects,nTrials];     // 1st choices, 1 or 2
  int<lower=1,upper=2> choice2[nSubjects,nTrials];     // 2nd choices, 1 or 2
  int<lower=0,upper=1> chswtch[nSubjects,nTrials];     // choice switch or not, 0 or 1
  real<lower=-1,upper=1> reward[nSubjects,nTrials];    // outcome, 1 or -1
  real<lower=0,upper=1>  with[nSubjects,nTrials];      // No. of with / 4
  real<lower=0,upper=1>  against[nSubjects,nTrials];   // No. of against / 4
}

transformed data {
  vector[2] initV;  // initial values for V
  initV <- rep_vector(0.0,2);    
}

parameters {
  // group-level parameters
  real myLR_mu_pr[2];
  real otherLR_mu_pr[2];
  real tau_mu_pr[2];
  real w_mu_pr[2];
  real coha_mu;
  real cohw_mu;
      
  real<lower=0> myLR_sd[2];
  real<lower=0> otherLR_sd[2];
  real<lower=0> tau_sd[2];
  real<lower=0> w_sd[2];
  real<lower=0> coha_sd;
  real<lower=0> cohw_sd;
    
  // subject-level raw parameters, follows norm(0,1), for later Matt Trick
  real myLR_raw[nSubjects,2];
  real otherLR_raw[nSubjects,2];
  real tau_raw[nSubjects,2];
  real w_raw[nSubjects,2];
  vector[nSubjects] coha_raw;
  vector[nSubjects] cohw_raw;
}

transformed parameters {
  // subject-level parameters
  real<lower=0,upper=1> myLR[nSubjects,2];
  real<lower=0,upper=1> otherLR[nSubjects,2];
  real<lower=0,upper=5> tau[nSubjects,2];
  real<lower=0,upper=1> w[nSubjects,2];
  vector[nSubjects] coha;
  vector[nSubjects] cohw;
  
  // Matt Trick. note that Phi_approx() takes 'real' rather than 'vector'
  for (s in 1:nSubjects) {
    for (i in 1:2) {
        myLR[s,i]    <- Phi_approx( myLR_mu_pr[i] + myLR_sd[i] * myLR_raw[s,i] );
        otherLR[s,i] <- Phi_approx( otherLR_mu_pr[i] + otherLR_sd[i] * otherLR_raw[s,i] );
        tau[s,i]     <- Phi_approx( tau_mu_pr[i] + tau_sd[i] * tau_raw[s,i] ) * 5;
        w[s,i]       <- Phi_approx( w_mu_pr[i] + w_sd[i] * w_raw[s,i] );
    }
  }
  coha <- coha_mu + coha_sd * coha_raw;
  cohw <- cohw_mu + cohw_sd * cohw_raw;
}

model {
  // values and PEs
  vector[2] choiceV[nTrials+1];    
  vector[2] myV[nTrials];
  vector[4] otherV[nTrials];
    
    
  // hyperparameters
  myLR_mu_pr     ~ normal(0,1);
  otherLR_mu_pr  ~ normal(0,1);
  tau_mu_pr      ~ normal(0,1);
  w_mu_pr        ~ normal(0,1);
  coha_mu        ~ normal(0,1);
  cohw_mu        ~ normal(0,1);
    
  myLR_sd    ~ cauchy(0,5);
  otherLR_sd ~ cauchy(0,5);
  tau_sd     ~ cauchy(0,5);
  w_sd       ~ cauchy(0,5);
  coha_sd    ~ cauchy(0,5);
  cohw_sd    ~ cauchy(0,5);
    
  // Matt Trick
  myLR_raw    ~ normal(0,1);
  otherLR_raw ~ normal(0,1);
  tau_raw     ~ normal(0,1);
  w_raw       ~ normal(0,1);
  coha_raw    ~ normal(0,1);
  cohw_raw    ~ normal(0,1);
    
  // subject loop and trial loop
  for (s in 1:nSubjects) {
    choiceV[1] <- rep_vector(0.0,2);
    myV[1]     <- rep_vector(0.0,2);
    otherV[1]  <- rep_vector(0.0,4);
    
    for (t in 1:nTrials) {
        choice1[s,t] ~ categorical_logit(tau[s,1] * choiceV[t]);
        
        
        
        
        
        
        
        
        
        
    } // trial
  }   // subject
    
    
    
    
    
    
    
    
    
    
    
}




generated quantities {
    
    
}












































