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
  real lr1_mu_pr;
  real lr2_mu_pr;
  real tau1_mu_pr;
  real tau2_mu_pr;
  real coha_mu;
  real cohw_mu;

  real<lower=0.01> lr1_sd;
  real<lower=0.01> lr2_sd;
  real<lower=0.01> tau1_sd;
  real<lower=0.01> tau2_sd;
  //real<lower=0.01, upper=3> coha_sd;
  //real<lower=0.01, upper=3> cohw_sd;
  real<lower=0.01> coha_sd;
  real<lower=0.01> cohw_sd;
  
  // subject-level raw parameters, follows norm(0,1), for later Matt Trick
  vector[nSubjects] lr1_raw;
  vector[nSubjects] lr2_raw;
  vector[nSubjects] tau1_raw;
  vector[nSubjects] tau2_raw;
  vector[nSubjects] coha_raw;
  vector[nSubjects] cohw_raw;
}

transformed parameters {
  // subject-level parameters
  vector<lower=0,upper=1>[nSubjects] lr1;
  vector<lower=0,upper=1>[nSubjects] lr2;
  vector<lower=0,upper=3>[nSubjects] tau1;
  vector<lower=0,upper=3>[nSubjects] tau2;
  vector[nSubjects] coha;
  vector[nSubjects] cohw;

  // Matt Trick. note that the input of Phi_approx must be 'real' rather than 'vector'
  for (s in 1:nSubjects) {
    lr1[s]  <- Phi_approx( lr1_mu_pr + lr1_sd * lr1_raw[s] );
    lr2[s]  <- Phi_approx( lr2_mu_pr + lr2_sd * lr2_raw[s] );
    tau1[s] <- Phi_approx( tau1_mu_pr + tau1_sd * tau1_raw[s] ) * 3;
    tau2[s] <- Phi_approx( tau2_mu_pr + tau2_sd * tau2_raw[s] ) * 3;
  }
  coha <- coha_mu + coha_sd * coha_raw;
  cohw <- cohw_mu + cohw_sd * cohw_raw;
}

model {
  // define the value and pe vectors
  vector[2] v1[nTrials+1];    // 1st values
  vector[2] v2[nTrials];      // 2nd values
  vector[nTrials] pe;         // prediction errors
  vector[nTrials] penc;       // pe for the non-chosen choice
  vector[2] tmp; 

  // hyperparameters
  lr1_mu_pr   ~ normal(0,1);
  lr2_mu_pr   ~ normal(0,1);
  tau1_mu_pr  ~ normal(0,1);
  tau2_mu_pr  ~ normal(0,1);
  coha_mu    ~ normal(0,1);
  cohw_mu    ~ normal(0,1);

  lr1_sd   ~ cauchy(0,5);
  lr2_sd   ~ cauchy(0,5);
  tau1_sd  ~ cauchy(0,5);
  tau2_sd  ~ cauchy(0,5);
  coha_sd ~ cauchy(0,3);
  cohw_sd ~ cauchy(0,3);
  
  // Matt Trick
  lr1_raw   ~ normal(0,1);
  lr2_raw   ~ normal(0,1);
  tau1_raw  ~ normal(0,1);
  tau2_raw  ~ normal(0,1);
  coha_raw ~ normal(0,1);
  cohw_raw ~ normal(0,1);

  // subject loop and trial loop
  for (s in 1:nSubjects) {
    v1[1] <- initV;
    
    for (t in 1:nTrials) {
      choice1[s,t] ~ categorical_logit(tau1[s] * v1[t]);
      
      v2[t][3-choice1[s,t]] <- v1[t][3-choice1[s,t]] + coha[s] * against[s,t];
      v2[t][choice1[s,t]]   <- v1[t][choice1[s,t]]   + cohw[s] * with[s,t];
    
      tmp <- tau2[s] * v2[t];
    
      //* compute action probs using built-in softmax function and related to choice data */
      choice2[s,t] ~ categorical_logit(tmp);
    
      //* prediction error */
      pe[t]   <-  reward[s,t] - v2[t][choice2[s,t]];
      penc[t] <- -reward[s,t] - v2[t][3-choice2[s,t]];

      //* value updating (learning) */
      if (choice1[s,t]==choice2[s,t]) {
          v1[t+1] <- v2[t] + lr1[s] * pe[t];
      } else {
          v1[t+1][choice2[s,t]]   <- v2[t][choice2[s,t]]   + lr1[s] * pe[t];
          v1[t+1][3-choice2[s,t]] <- v2[t][3-choice2[s,t]] + lr2[s] * penc[t];  
      } // if
    }
  }
}

generated quantities {
  real<lower=0,upper=1> lr1_mu;
  real<lower=0,upper=1> lr2_mu;
  real<lower=0,upper=3> tau1_mu;
  real<lower=0,upper=3> tau2_mu;
  
  real log_lik1[nSubjects]; 
  real log_lik2[nSubjects]; 
  vector[2] ev1[nTrials+1];
  vector[2] ev2[nTrials];
  vector[nTrials] pe2;
  vector[nTrials] penc2;
  
  lr1_mu  <- Phi_approx(lr1_mu_pr);
  lr2_mu  <- Phi_approx(lr2_mu_pr);
  tau1_mu <- Phi_approx(tau1_mu_pr) * 3;
  tau2_mu <- Phi_approx(tau2_mu_pr) * 3;
  
  for (s in 1:nSubjects) {
    log_lik1[s] <- 0;
    log_lik2[s] <- 0;
    ev1[1] <- rep_vector(0.0,2);
    
    for (t in 1:nTrials) {
      log_lik1[s] <- log_lik1[s] + categorical_logit_log(choice1[s,t], tau1[s] * ev1[t]);

      ev2[t][3-choice1[s,t]] <- ev1[t][3-choice1[s,t]] + coha[s] * against[s,t];
      ev2[t][choice1[s,t]]   <- ev1[t][choice1[s,t]]   + cohw[s] * with[s,t];
      
      log_lik2[s] <- log_lik2[s] + categorical_logit_log(choice2[s,t], tau2[s] * ev2[t]);
      
      pe2[t]   <-  reward[s,t] - ev2[t][choice2[s,t]];
      penc2[t] <- -reward[s,t] - ev2[t][3-choice2[s,t]];

      if (choice1[s,t]==choice2[s,t]) {
          ev1[t+1] <- ev2[t] + lr1[s] * pe2[t];
      } else {
          ev1[t+1][choice2[s,t]]   <- ev2[t][choice2[s,t]]   + lr1[s] * pe2[t];
          ev1[t+1][3-choice2[s,t]] <- ev2[t][3-choice2[s,t]] + lr2[s] * penc2[t];  
      }
    }
  }
}