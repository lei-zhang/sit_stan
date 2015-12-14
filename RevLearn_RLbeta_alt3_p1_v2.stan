data {
  int<lower=1> nSubjects;                                    // number of subjects
  int<lower=1> nTrials;                                      // number of trials 
  int<lower=1,upper=2> otherChoice2[nSubjects,nTrials,4];
  real<lower=-1,upper=1> otherReward[nSubjects,nTrials,4];
}

transformed data {
  vector[2] initV; 
  initV <- rep_vector(0.0,2);  
}

parameters {
  // group-level parameters
  vector[4] lr_mu_pr;
  vector[4] tau_mu_pr; 
  vector<lower=0>[4] lr_sd;
  vector<lower=0>[4] tau_sd;
  
  // subject-level raw parameters, follows norm(0,1), for later Matt Trick
  vector[nSubjects] lr_raw[4];        // dim: [4 nSubjects]
  vector[nSubjects] tau_raw[4];        // dim: [4 nSubjects]  
}

transformed parameters {
  // subject-level parameters
  vector<lower=0,upper=1>[nSubjects] lr[4];
  vector<lower=0,upper=5>[nSubjects] tau[4];
  
  // Matt Trick, note that the input of Phi_approx must be 'real' rather than 'vector'
  for (o in 1:4)           {
    for (s in 1:nSubjects) {
      lr[o,s]    <- Phi_approx( lr_mu_pr[o] + lr_sd[o] * lr_raw[o,s] );
      tau[o,s]   <- Phi_approx( tau_mu_pr[o] + tau_sd[o] * tau_raw[o,s] ) * 5;      
    }
  }
}

model {
  // define the value and pe vectors
  vector[2] v[nTrials+1,4];  // dim [101 4 2]
  real pe[nTrials,4];
  
  // hyperparameters
  lr_mu_pr  ~ normal(0,1);
  tau_mu_pr ~ normal(0,1);
  lr_sd     ~ cauchy(0,5);
  tau_sd    ~ cauchy(0,5);
  
  // Matt Trick
  for (o in 1:4) {
    lr_raw[o]  ~ normal(0,1);
    tau_raw[o] ~ normal(0,1);
  }

  // subject loop and trial loop
  for (s in 1:nSubjects) {
    for (o in 1:4) {
      v[1,o] <- initV;
    }
          
    for (t in 1:nTrials) {
      for (o in 1:4) {
         otherChoice2[s,t,o] ~ categorical_logit( tau[o,s] * v[t,o] );
         pe[t,o]  <- otherReward[s,t,o] - v[t,o,otherChoice2[s,t,o]];
         v[t+1,o] <- v[t,o];
         v[t+1,o,otherChoice2[s,t,o]] <- v[t,o,otherChoice2[s,t,o]] + lr[o,s] * pe[t,o]; 
      }
    }  // trial loop
  }    // subject loop
}

generated quantities {
  real<lower=0,upper=1> lr_mu[4];
  real<lower=0,upper=5> tau_mu[4];
  
  for (o in 1:4) {
    lr_mu[o]  <- Phi_approx(lr_mu_pr[o]);
    tau_mu[o] <- Phi_approx(tau_mu_pr[o])*5;
  }
}