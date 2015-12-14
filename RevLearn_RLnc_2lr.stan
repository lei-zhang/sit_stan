data {
  int<lower=1> nSubjects;                              // number of subjects
  int<lower=1> nTrials;                                // number of trials 
  int<lower=1,upper=2> choice2[nSubjects,nTrials];      // 2nd choices, 1 or 2
  real<lower=-1,upper=1> reward[nSubjects,nTrials];    // outcome, 1 or -1
}

transformed data {
  vector[2] initV;  // initial values for V
  initV <- rep_vector(0.0,2);    
}

parameters {
  // group-level parameters
  real lr1_mu_pr;    // lr1_mu before probit
  real lr2_mu_pr;    // lr2_mu before probit
  real tau_mu_pr;    // tau_mu before probit  

  real<lower=0> lr1_sd;
  real<lower=0> lr2_sd;
  real<lower=0> tau_sd;
  
  // subject-level raw parameters, follows norm(0,1), for later Matt Trick
  vector[nSubjects] lr1_raw;
  vector[nSubjects] lr2_raw;
  vector[nSubjects] tau_raw;
}

transformed parameters {
  // subject-level parameters
  vector<lower=0,upper=1>[nSubjects] lr1;
  vector<lower=0,upper=1>[nSubjects] lr2;
  vector<lower=0,upper=10>[nSubjects] tau;
  
  // Matt Trick, note that the input of Phi_approx must be 'real' rather than 'vector'
  for (s in 1:nSubjects) {
    lr1[s] <- Phi_approx( lr1_mu_pr + lr1_sd * lr1_raw[s] );
    lr2[s] <- Phi_approx( lr2_mu_pr + lr2_sd * lr2_raw[s] );
    tau[s] <- Phi_approx( tau_mu_pr + tau_sd * tau_raw[s] ) * 10;
  }
}

model {
  // define the value and pe vectors
  vector[2] v[nTrials+1];     // values
  vector[nTrials] pe;         // prediction errors
  vector[nTrials] penc;       // pe for the non-chosen choice

  // hyperparameters
  lr1_mu_pr ~ normal(0,1);
  lr2_mu_pr ~ normal(0,1);
  tau_mu_pr ~ normal(0,1);
  
  lr1_sd ~ cauchy(0,5);
  lr2_sd ~ cauchy(0,5);
  tau_sd ~ cauchy(0,5);
  
  // Matt Trick
  lr1_raw ~ normal(0,1);
  lr2_raw ~ normal(0,1);
  tau_raw ~ normal(0,1);
  
  // subject loop and trial loop
  for (s in 1:nSubjects) {
    v[1] <- initV;  

    for (t in 1:nTrials) {
            
      //* compute action probs using built-in softmax function and related to choice data */
      choice2[s,t] ~ categorical_logit( tau[s] * v[t] );

      //* prediction error */
      pe[t]   <-  reward[s,t] - v[t][choice2[s,t]];
      penc[t] <- -reward[s,t] - v[t][3-choice2[s,t]];

      //* value updating (learning) */
      v[t+1][choice2[s,t]]   <- v[t][choice2[s,t]]   + lr1[s] * pe[t];
      v[t+1][3-choice2[s,t]] <- v[t][3-choice2[s,t]] + lr2[s] * penc[t];
    }
  }
}

generated quantities {
  real<lower=0,upper=1> lr1_mu;
  real<lower=0,upper=1> lr2_mu;
  real<lower=0,upper=10> tau_mu;
  
  real log_lik[nSubjects]; 
  vector[2] v2[nTrials+1];
  vector[nTrials] pe2;
  vector[nTrials] penc2;
  int<lower=1,upper=2> c_rep[nSubjects, nTrials];
  
  lr1_mu <- Phi_approx(lr1_mu_pr);
  lr2_mu <- Phi_approx(lr2_mu_pr);
  tau_mu <- Phi_approx(tau_mu_pr) * 10;
  
  for (s in 1:nSubjects) {
    log_lik[s] <- 0;
    v2[1] <- initV;  

    for (t in 1:nTrials) {
      log_lik[s] <- log_lik[s] + categorical_logit_log(choice2[s,t], tau[s] * v2[t]);
      c_rep[s,t] <- categorical_rng( softmax(tau[s]*v2[t]) );

      pe2[t]   <-  reward[s,t] - v2[t][choice2[s,t]];
      penc2[t] <- -reward[s,t] - v2[t][3-choice2[s,t]];
      v2[t+1][choice2[s,t]]   <- v2[t][choice2[s,t]]   + lr1[s] * pe2[t];
      v2[t+1][3-choice2[s,t]] <- v2[t][3-choice2[s,t]] + lr2[s] * penc2[t];
    }
  }
}