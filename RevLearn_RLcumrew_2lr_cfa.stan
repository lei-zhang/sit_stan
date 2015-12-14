data {
  int<lower=1> nSubjects;                              // number of subjects
  int<lower=1> nTrials;                                // number of trials 
  int<lower=1,upper=2> choice1[nSubjects,nTrials];     // 1st choices, 1 or 2
  int<lower=1,upper=2> choice2[nSubjects,nTrials];     // 2nd choices, 1 or 2
  real<lower=-1,upper=1> reward[nSubjects,nTrials];    // outcome, 1 or -1
  real<lower=0,upper=1>  otherReward[nSubjects,nTrials,4];
  real<lower=0,upper=1>  otherWith[nSubjects,nTrials,4];
}

transformed data {
  vector[2] initV;       // initial values for V
  row_vector[4] initCR;  // initial values for cr
  
  initV  <- rep_vector(0.0,2);    
  initCR <- rep_row_vector(0.25,4);
}

parameters {
  // group-level parameters
  real lr1_mu_pr;   // lr1_mu before probit
  real lr2_mu_pr;   // lr2_mu before probit
  real tau_mu_pr;   // tau_mu before probit  
  real disc_mu_pr;  // discounting gamma, before probit
  real cfa_mu_pr;
  real cra_mu;
  real crw_mu;

  real<lower=0> lr1_sd;
  real<lower=0> lr2_sd;
  real<lower=0> tau_sd;
  real<lower=0> disc_sd;
  real<lower=0> cfa_sd;
  real<lower=0> cra_sd;
  real<lower=0> crw_sd;
  
  // subject-level row parameters, follows norm(0,1), for later Matt Trick
  vector[nSubjects] lr1_raw;
  vector[nSubjects] lr2_raw;
  vector[nSubjects] tau_raw;
  vector[nSubjects] disc_raw;
  vector[nSubjects] cfa_raw;
  vector[nSubjects] cra_raw;
  vector[nSubjects] crw_raw;
}

transformed parameters {
  // subject-level parameters
  vector<lower=0,upper=1>[nSubjects] lr1;
  vector<lower=0,upper=1>[nSubjects] lr2;
  vector<lower=0,upper=10>[nSubjects] tau;
  vector<lower=0,upper=1>[nSubjects]  disc;
  vector<lower=0,upper=1>[nSubjects]  cfa;
  vector[nSubjects] cra;
  vector[nSubjects] crw;

  // Matt Trick, note that the input of Phi_approx must be 'real' rather than 'vector'
  for (s in 1:nSubjects) {
    lr1[s] <- Phi_approx( lr1_mu_pr + lr1_sd * lr1_raw[s] );
    lr2[s] <- Phi_approx( lr2_mu_pr + lr2_sd * lr2_raw[s] );
    tau[s]  <- Phi_approx( tau_mu_pr  + tau_sd  * tau_raw[s] ) * 10;
    disc[s] <- Phi_approx( disc_mu_pr + disc_sd * disc_raw[s] );
    cfa[s]  <- Phi_approx( cfa_mu_pr + cfa_sd * cfa_raw[s] );
  }
  cra <- cra_mu + cra_sd * cra_raw; // vectorization
  crw <- crw_mu + crw_sd * crw_raw;
}

model {
  // define the value and pe vectors
  vector[2] v[nTrials+1];  // values
  vector[nTrials] pe;      // prediction errors
  vector[nTrials] penc;    // pe for the non-chosen choice  
  matrix[nTrials,4] cr;    // cumulative reward
  
  // hyperparameters
  lr1_mu_pr  ~ normal(0,1);
  lr2_mu_pr  ~ normal(0,1);
  tau_mu_pr  ~ normal(0,1);
  disc_mu_pr ~ normal(0,1);
  cfa_mu_pr  ~ normal(0,1);
  cra_mu     ~ normal(0,1);
  crw_mu     ~ normal(0,1);
  
  lr1_sd  ~ cauchy(0,5);
  lr2_sd  ~ cauchy(0,5);
  tau_sd  ~ cauchy(0,5);
  disc_sd ~ cauchy(0,5);
  cra_sd  ~ cauchy(0,5);
  crw_sd  ~ cauchy(0,5);
  cfa_sd  ~ cauchy(0,5);
  
  // Matt Trick
  lr1_raw  ~ normal(0,1);
  lr2_raw  ~ normal(0,1);
  tau_raw  ~ normal(0,1);
  disc_raw ~ normal(0,1);
  cfa_raw  ~ normal(0,1);
  cra_raw  ~ normal(0,1);
  crw_raw  ~ normal(0,1);
  
  // subject loop and trial loop
  for (s in 1:nSubjects) {
    v[1] <- initV;
    
    for (t in 1:nTrials) {      
      // start calculating cumulative reward --------------------
      matrix[t-1,4] wgtrew;
      
      if (t <= 2) {
        cr[t] <- initCR;
      } else {
        for (ct in 1:(t-1)) {
          wgtrew[ct] <- (disc[s] ^ (t-ct)) * to_row_vector(otherReward[s,ct]);
        }
        for (o in 1:4) {
          cr[t,o] <- sum( sub_col(wgtrew, 1,o,t-1));
        }
      }
      // cr calculation finished -------------------------

      // re-weight value after have seen the group decisions
      v[t][3-choice1[s,t]] <- v[t][3-choice1[s,t]] +cra[s] *sum( cr[t]/(sum(cr[t])+machine_precision()) .* (-1*to_row_vector(otherWith[s,t])+1) );
      v[t][choice1[s,t]]   <- v[t][choice1[s,t]]   +crw[s] *sum( cr[t]/(sum(cr[t])+machine_precision()) .* to_row_vector(otherWith[s,t]) );

      //* compute action probs using built-in softmax function and related to choice data */
      choice2[s,t] ~ categorical_logit( tau[s] * v[t] );

      //* prediction error */
      pe[t]   <-  reward[s,t] - v[t][choice2[s,t]];
      penc[t] <- (-reward[s,t]*cfa[s]) - v[t][3-choice2[s,t]];

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
  real<lower=0,upper=1>  disc_mu;
  real<lower=0,upper=1>  cfa_mu;

  real log_lik[nSubjects]; 
  vector[2] v2[nTrials+1];
  vector[nTrials] pe2;
  vector[nTrials] penc2;
  matrix[nTrials,4] cr2;
  int<lower=1,upper=2> c_rep[nSubjects, nTrials];
  
  lr1_mu  <- Phi_approx(lr1_mu_pr);
  lr2_mu  <- Phi_approx(lr2_mu_pr);
  tau_mu  <- Phi_approx(tau_mu_pr) * 10;
  disc_mu <- Phi_approx( disc_mu_pr );
  cfa_mu  <- Phi_approx(cfa_mu_pr);

  for (s in 1:nSubjects) {
    log_lik[s] <- 0;
    v2[1] <- initV;
    
    for (t in 1:nTrials) {      
      matrix[t-1,4] wgtrew2;
      
      if (t <= 2) {
        cr2[t] <- initCR;
      } else {
        for (ct in 1:(t-1)) {
          wgtrew2[ct] <- (disc[s] ^ (t-ct)) * to_row_vector(otherReward[s,ct]);
        }
        for (o in 1:4) {
          cr2[t,o] <- sum( sub_col(wgtrew2, 1,o,t-1));
        }
      }

      v2[t][3-choice1[s,t]] <- v2[t][3-choice1[s,t]] + cra[s] * sum( cr2[t]/(sum(cr2[t])+machine_precision()) .* (-1*to_row_vector(otherWith[s,t])+1) );
      v2[t][choice1[s,t]]   <- v2[t][choice1[s,t]]   + crw[s] * sum( cr2[t]/(sum(cr2[t])+machine_precision()) .* to_row_vector(otherWith[s,t]) );

      log_lik[s] <- log_lik[s] + categorical_logit_log(choice2[s,t], tau[s] * v2[t]);
      c_rep[s,t] <- categorical_rng( softmax(tau[s]*v2[t]) );

      pe2[t]   <-  reward[s,t] - v2[t][choice2[s,t]];
      penc2[t] <- (-reward[s,t]*cfa[s]) - v2[t][3-choice2[s,t]];

      v2[t+1][choice2[s,t]]   <- v2[t][choice2[s,t]]   + lr1[s] * pe2[t];
      v2[t+1][3-choice2[s,t]] <- v2[t][3-choice2[s,t]] + lr2[s] * penc2[t];
    }
  }
}