data {
  int<lower=1> nSubjects;                              // number of subjects
  int<lower=1> nTrials;                                // number of trials 
  int<lower=1,upper=2> choice2[nSubjects,nTrials];     // 2nd choices, 1 or 2, must be 'int'
  real<lower=-1,upper=1> reward[nSubjects,nTrials];    // outcome, 1 or -1, 'real' is faster than 'int'
}

transformed data {
  vector[2] initV;  // initial values for V
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
  vector<lower=0,upper=10>[nSubjects] tau;
  
  // Matt Trick, note that the input of Phi_approx must be 'real' rather than 'vector'
  for (s in 1:nSubjects) {
    lr[s]  <- Phi_approx( lr_mu_pr + lr_sd * lr_raw[s] );
    tau[s] <- Phi_approx( tau_mu_pr + tau_sd * tau_raw[s] )*10;
  }
}

model {
  // define the value and pe vectors
  vector[2] v[nTrials+1];     // values
  vector[nTrials] pe;         // prediction errors  

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
    v[1] <- initV;  

    for (t in 1:nTrials) {
      //* compute action probs using built-in softmax function and related to choice data */
      // using chpice_minus_one and apply bernoulli_logit() makes the sampling even longer...
      choice2[s,t] ~ categorical_logit( tau[s] * v[t] );

      //* prediction error */
      pe[t] <- reward[s,t] - v[t][choice2[s,t]];

      //* value updating (learning) */
      v[t+1] <- v[t]; // make a copy of current value into t+1
      v[t+1][choice2[s,t]] <- v[t][choice2[s,t]] + lr[s] * pe[t]; // overwrite chosen value with pe update
    }
  }
}

generated quantities {
  real<lower=0,upper=1> lr_mu; 
  real<lower=0,upper=10> tau_mu;
  
  real log_lik[nSubjects]; 
  vector[2] v2[nTrials+1];
  vector[nTrials] pe2;
  int<lower=1,upper=2> c_rep[nSubjects, nTrials];

  // recover the mu and sigma
  lr_mu   <- Phi_approx(lr_mu_pr);
  tau_mu  <- Phi_approx(tau_mu_pr)*10;

  // compute the log-likelihood and reproduce data for posterior predictive check
  for (s in 1:nSubjects) {
    v2[1] <- initV;
    log_lik[s] <- 0;
   
    for (t in 1:nTrials) {
      log_lik[s] <- log_lik[s] + categorical_logit_log(choice2[s,t], tau[s] * v2[t]);
      c_rep[s,t] <- categorical_rng( softmax(tau[s]*v2[t]) );
      
      pe2[t]  <- reward[s,t] - v2[t][choice2[s,t]];
      v2[t+1] <- v2[t];
      v2[t+1][choice2[s,t]] <- v2[t][choice2[s,t]] + lr[s] * pe2[t];
    }
  }
}