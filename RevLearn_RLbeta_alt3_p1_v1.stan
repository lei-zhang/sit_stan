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
  real lr_mu_pr;    // lr_mu before probit
  real tau_mu_pr;   // tau_mu before probit  
  real<lower=0> lr_sd;
  real<lower=0> tau_sd;
  
  // subject-level raw parameters, follows norm(0,1), for later Matt Trick
  vector[nSubjects] lr_raw;
  vector[nSubjects] tau_raw;
}

transformed parameters {
  // subject-level parameters
  vector<lower=0,upper=1>[nSubjects] lr;
  vector<lower=0,upper=5>[nSubjects] tau;

  // Matt Trick, note that the input of Phi_approx must be 'real' rather than 'vector'
  for (s in 1:nSubjects) {
    lr[s]  <- Phi_approx( lr_mu_pr + lr_sd * lr_raw[s] );
    tau[s] <- Phi_approx( tau_mu_pr + tau_sd * tau_raw[s] ) * 5;
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
  lr_raw  ~ normal(0,1);
  tau_raw ~ normal(0,1);
  
  // subject loop and trial loop
  for (s in 1:nSubjects) {
    for (o in 1:4) {
      v[1,o] <- initV;
    }
          
    for (t in 1:nTrials) {
      for (o in 1:4) {
         otherChoice2[s,t,o] ~ categorical_logit( tau[s] * v[t,o] );
         pe[t,o]  <- otherReward[s,t,o] - v[t,o,otherChoice2[s,t,o]];
         v[t+1,o] <- v[t,o];
         v[t+1,o,otherChoice2[s,t,o]] <- v[t,o,otherChoice2[s,t,o]] + lr[s] * pe[t,o]; 
      }
    }  // trial loop
  }    // subject loop
}

generated quantities {
  real<lower=0,upper=1> lr_mu; 
  real<lower=0,upper=5> tau_mu;

  lr_mu   <- Phi_approx(lr_mu_pr);
  tau_mu  <- Phi_approx(tau_mu_pr)*5;
}