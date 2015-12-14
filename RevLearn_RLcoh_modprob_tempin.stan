data {
  int<lower=1> nSubjects;                              // number of subjects
  int<lower=1> nTrials;                                // number of trials 
  int<lower=1,upper=2> choice1[nSubjects,nTrials];     // 1st choices, 1 or 2
  int<lower=1,upper=2> choice2[nSubjects,nTrials];     // 2nd choices, 1 or 2
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
  real lr_mu_pr;    // lr_mu before probit
  real tau_mu_pr;   // tau_mu before probit  
  real coha_mu;
  real cohw_mu;

  real<lower=0.01> lr_sd;
  real<lower=0.01> tau_sd;
  //real<lower=0.01, upper=3> coha_sd;
  //real<lower=0.01, upper=3> cohw_sd;
  real<lower=0.01> coha_sd;
  real<lower=0.01> cohw_sd;
  
  // subject-level raw parameters, follows norm(0,1), for later Matt Trick
  vector[nSubjects] lr_raw;
  vector[nSubjects] tau_raw;
  vector[nSubjects] coha_raw;
  vector[nSubjects] cohw_raw;
}

transformed parameters {
  // subject-level parameters
  vector<lower=0,upper=1>[nSubjects] lr;
  vector<lower=0,upper=2.5>[nSubjects] tau;
  vector[nSubjects] coha;
  vector[nSubjects] cohw;

  // Matt Trick. note that the input of Phi_approx must be 'real' rather than 'vector'
  for (s in 1:nSubjects) {
    lr[s]  <- Phi_approx( lr_mu_pr + lr_sd * lr_raw[s] );
    tau[s] <- Phi_approx( tau_mu_pr + tau_sd * tau_raw[s] ) * 2.5;
  }
  coha <- coha_mu + coha_sd * coha_raw;
  cohw <- cohw_mu + cohw_sd * cohw_raw;
}

model {
  // define the value and pe vectors
  vector[2] v1[nTrials+1];    // 1st values
  vector[2] tmpv;             // temp values
  vector[nTrials] pe;         // prediction errors
  vector[nTrials] penc;       // pe for the non-chosen choice

  // hyperparameters
  lr_mu_pr   ~ normal(0,1);
  tau_mu_pr  ~ normal(0,1);
  coha_mu    ~ normal(0,1);
  cohw_mu    ~ normal(0,1);

  lr_sd   ~ cauchy(0,5);
  tau_sd  ~ cauchy(0,5);
  coha_sd ~ cauchy(0,3);
  cohw_sd ~ cauchy(0,3);
  
  // Matt Trick
  lr_raw   ~ normal(0,1);
  tau_raw  ~ normal(0,1);
  coha_raw ~ normal(0,1);
  cohw_raw ~ normal(0,1);

  // subject loop and trial loop
  for (s in 1:nSubjects) {
    v1[1] <- initV;
    
    for (t in 1:nTrials) {
      choice1[s,t] ~ categorical_logit(tau[s] * v1[t]);
      
      tmpv[3-choice1[s,t]] <- tau[s] * v1[t][3-choice1[s,t]] + coha[s] * against[s,t];
      tmpv[choice1[s,t]]   <- tau[s] * v1[t][choice1[s,t]]   + cohw[s] * with[s,t];
    
      //* compute action probs using built-in softmax function and related to choice data */
      choice2[s,t] ~ categorical_logit(tmpv);
    
      //* prediction error */
      pe[t]   <-  reward[s,t] - v1[t][choice2[s,t]];
      penc[t] <- -reward[s,t] - v1[t][3-choice2[s,t]];

      //* value updating (learning) */
      if (choice1[s,t]==choice2[s,t]) {
          v1[t+1] <- v1[t] + lr[s] * pe[t];
      } else {
          v1[t+1][choice2[s,t]]   <- v1[t][choice2[s,t]]   + lr[s] * pe[t];
          v1[t+1][3-choice2[s,t]] <- v1[t][3-choice2[s,t]] + lr[s] * penc[t];  
      } // if
    }
  }
}

generated quantities {
  real<lower=0,upper=1> lr_mu; 
  real<lower=0,upper=2.5> tau_mu;
  
  real log_lik1[nSubjects]; 
  real log_lik2[nSubjects]; 
  vector[2] ev1[nTrials+1];
  vector[2] tmpv1;
  vector[nTrials] pe2;
  vector[nTrials] penc2;
  
  lr_mu  <- Phi_approx(lr_mu_pr);
  tau_mu <- Phi_approx(tau_mu_pr) * 2.5;
  
  for (s in 1:nSubjects) {
    log_lik1[s] <- 0;
    log_lik2[s] <- 0;
    ev1[1] <- rep_vector(0.0,2);
    
    for (t in 1:nTrials) {
      log_lik1[s] <- log_lik1[s] + categorical_logit_log(choice1[s,t], tau[s] * ev1[t]);
      
      tmpv1[3-choice1[s,t]] <- tau[s] * ev1[t][3-choice1[s,t]] + coha[s] * against[s,t];
      tmpv1[choice1[s,t]]   <- tau[s] * ev1[t][choice1[s,t]]   + cohw[s] * with[s,t];
      
      log_lik2[s] <- log_lik2[s] + categorical_logit_log(choice2[s,t], tmpv1);
      
      pe2[t]   <-  reward[s,t] - ev1[t][choice2[s,t]];
      penc2[t] <- -reward[s,t] - ev1[t][3-choice2[s,t]];

      if (choice1[s,t]==choice2[s,t]) {
          ev1[t+1] <- ev1[t] + lr[s] * pe2[t];
      } else {
          ev1[t+1][choice2[s,t]]   <- ev1[t][choice2[s,t]]   + lr[s] * pe2[t];
          ev1[t+1][3-choice2[s,t]] <- ev1[t][3-choice2[s,t]] + lr[s] * penc2[t];  
      }
    }
  }
}