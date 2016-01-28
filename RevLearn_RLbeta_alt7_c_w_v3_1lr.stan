data {
  int<lower=1> nSubjects;                                  // number of subjects
  int<lower=1> nTrials;                                    // number of trials 
  int<lower=1,upper=2> choice1[nSubjects,nTrials];         // 1st choices, 1 or 2
  int<lower=1,upper=2> choice2[nSubjects,nTrials];         // 2nd choices, 1 or 2
  int<lower=0,upper=1> chswtch[nSubjects,nTrials];         // choice switch, 0 or 1
  real<lower=-1,upper=1> reward[nSubjects,nTrials];        // outcome, 1 or -1
  real<lower=0,upper=1>  wgtWith[nSubjects,nTrials];       
  real<lower=0,upper=1>  wgtAgst[nSubjects,nTrials];   
  matrix<lower=0,upper=1>[nTrials,4]  otherReward2[nSubjects]; // dim: [ns, t, 4]
  matrix<lower=-1,upper=1>[nTrials,4] otherReward2_actural[nSubjects]; 
  matrix<lower=0,upper=1>[nTrials,4]  otherWith2[nSubjects];
  matrix<lower=0.25,upper=1>[nTrials,4]  wOthers[nSubjects]; 
}

transformed data {
  vector[2] initV;    // initial values for V
  int B;              // number of predictors
  vector[3] pwr;      // power
  
  initV  <- rep_vector(0.0,2); 
  B      <- 4;
  pwr[1] <- 2;
  pwr[2] <- 1;
  pwr[3] <- 0;
}

parameters {
  // group-level parameters
  real lr_mu_pr;
  real disc_mu_pr; 
  real beta0_mu_pr;
  vector[B] beta_mu;

  real<lower=0> lr_sd;
  real<lower=0> disc_sd;
  real<lower=0> beta0_sd;
  vector<lower=0>[B] beta_sd;
  
  // subject-level raw parameters, follows norm(0,1), for later Matt Trick
  vector[nSubjects] lr_raw;        // dim: [1, nSubjects]
  vector[nSubjects] disc_raw;
  vector[nSubjects] beta0_raw;
  vector[nSubjects] beta_raw[B];   // dim: [B, nSubjects]
}

transformed parameters {
  // subject-level parameters
  vector<lower=0,upper=1>[nSubjects] lr;
  vector<lower=0,upper=1+1e-10>[nSubjects] disc;
  vector<lower=0,upper=2>[nSubjects] beta0;
  vector[nSubjects] beta[B];
  
  // Matt Trick, note that the input of Phi_approx must be 'real' rather than 'vector'
  for (s in 1:nSubjects) {
    lr[s]    <- Phi_approx( lr_mu_pr + lr_sd * lr_raw[s] );
    disc[s]  <- Phi_approx( disc_mu_pr + disc_sd * disc_raw[s] ) + machine_precision();
    beta0[s] <- Phi_approx( beta0_mu_pr + beta0_sd * beta0_raw[s] ) * 2;
  }
  for (i in 1:B) {
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
  matrix[3,4] disc_mat;
  matrix[3,4] reward_mat;
  row_vector[4] othW;
  row_vector[4] otherReward_vec;
  
  // hyperparameters
  lr_mu_pr   ~ normal(0,1);
  disc_mu_pr ~ normal(0,1);
  beta_mu    ~ normal(0,1);
  
  lr_sd   ~ cauchy(0,5);
  disc_sd ~ cauchy(0,5);
  beta_sd ~ cauchy(0,5);
  
  // Matt Trick
  lr_raw   ~ normal(0,1);
  disc_raw ~ normal(0,1);
  
  for (i in 1:B) {
    beta_raw[i] ~ normal(0,1);
  } 
  
  // subject loop and trial loop
  for (s in 1:nSubjects) {
    myValue[1]    <- initV;
    otherValue[1] <- initV;

    for (t in 1:nTrials) {
      valfun1 <- beta0[s] * myValue[t] + (2-beta0[s]) * otherValue[t];
      choice1[s,t] ~ categorical_logit( valfun1 );

      valdiff <- valfun1[choice1[s,t]] - valfun1[3-choice1[s,t]];
      valfun2 <- beta[1,s] + beta[2,s] * valdiff + beta[3,s] * wgtWith[s,t] + beta[4,s] * wgtAgst[s,t];
      chswtch[s,t] ~ bernoulli_logit(valfun2);

      // my prediction error
      pe[t]   <-  reward[s,t] - myValue[t,choice2[s,t]];
      penc[t] <- -reward[s,t] - myValue[t,3-choice2[s,t]];

      // update my value
      myValue[t+1,choice2[s,t]]   <- myValue[t,choice2[s,t]]   + lr[s] * pe[t];
      myValue[t+1,3-choice2[s,t]] <- myValue[t,3-choice2[s,t]] + lr[s] * penc[t];

      // use weighted discounted outcome to update the others' value
      othW <- otherWith2[s,t]; 

      if (t==1) {
        otherReward_vec  <- wOthers[s,t] .* otherReward2_actural[s,t];
      } else if (t==2) {
        otherReward_vec  <- wOthers[s,t] .* otherReward2_actural[s,t] + wOthers[s,t-1] .* otherReward2[s,t-1] * disc[s];
      } else {
        disc_mat <- rep_matrix(exp(log(disc[s])*pwr),4); // replicate by 4 columns 
        reward_mat <- append_row( block(otherReward2[s], t-2, 1, 2, 4), otherReward2_actural[s,t] );
        otherReward_vec  <- rep_row_vector(1.0, 3) * (disc_mat .* block(wOthers[s], t-2, 1, 3, 4) .* reward_mat);
      }

      otherValue[t+1,choice2[s,t]]   <- sum( othW     .* otherReward_vec );
      otherValue[t+1,3-choice2[s,t]] <- sum( (1-othW) .* otherReward_vec );
      otherValue[t+1] <- (1 ./ ( 1 + exp(-otherValue[t+1]))) * 2 - 1;

    }  // trial loop
  }    // subject loop
}

generated quantities {
  real<lower=0,upper=1> lr_mu;
  real<lower=0,upper=1+1e-10> disc_mu; 
  real<lower=0,upper=2> beta0_mu;
  
  real log_likc1[nSubjects]; 
  real log_likc2[nSubjects]; 
  vector[2] myValue2[nTrials+1];
  vector[2] otherValue2[nTrials+1];
  vector[nTrials] pe2;
  vector[nTrials] penc2;
  real valdiff_gen;
  vector[2] valfun1_gen;
  real valfun2_gen;
  matrix[3,4] disc_mat2;
  matrix[3,4] reward_mat2;
  row_vector[4] othW2;
  row_vector[4] otherReward_vec2;

  int<lower=1,upper=2> c1_rep[nSubjects, nTrials];
  int<lower=0,upper=1> c2_rep[nSubjects, nTrials];

  lr_mu    <- Phi_approx(lr_mu_pr);
  disc_mu  <- Phi_approx(disc_mu_pr) + machine_precision();
  beta0_mu <- Phi_approx(beta0_mu_pr) * 2;

  for (s in 1:nSubjects) {
    myValue2[1]    <- initV;
    otherValue2[1] <- initV;
    log_likc1[s]   <- 0;
    log_likc2[s]   <- 0;

    for (t in 1:nTrials) {
      valfun1_gen  <- beta0[s] * myValue2[t] + (2-beta0[s]) * otherValue2[t];
      log_likc1[s] <- log_likc1[s] + categorical_logit_log(choice1[s,t], valfun1_gen);

      c1_rep[s,t] <- categorical_rng( softmax(valfun1_gen) );

      valdiff_gen  <- valfun1_gen[choice1[s,t]] - valfun1_gen[3-choice1[s,t]];
      valfun2_gen  <- beta[1,s] + beta[2,s] * valdiff_gen + beta[3,s] * wgtWith[s,t] + beta[4,s] * wgtAgst[s,t];
      log_likc2[s] <- log_likc2[s] + bernoulli_logit_log(chswtch[s,t], valfun2_gen);
   
      c2_rep[s,t]   <- bernoulli_rng( inv_logit(valfun2_gen) );

      pe2[t]   <-  reward[s,t] - myValue2[t,choice2[s,t]];
      penc2[t] <- -reward[s,t] - myValue2[t,3-choice2[s,t]];
      
      myValue2[t+1,choice2[s,t]]   <- myValue2[t,choice2[s,t]]   + lr[s] * pe2[t];
      myValue2[t+1,3-choice2[s,t]] <- myValue2[t,3-choice2[s,t]] + lr[s] * penc2[t];

      othW2 <- otherWith2[s,t]; 

      if (t==1) {
        otherReward_vec2  <- wOthers[s,t] .* otherReward2_actural[s,t];
      } else if (t==2) {
        otherReward_vec2  <- wOthers[s,t] .* otherReward2_actural[s,t] + wOthers[s,t-1] .* otherReward2[s,t-1] * disc[s];
      } else {
        disc_mat2   <- rep_matrix(exp(log(disc[s])*pwr),4); // replicate by 4 columns 
        reward_mat2 <- append_row( block(otherReward2[s], t-2, 1, 2, 4), otherReward2_actural[s,t] );
        otherReward_vec2  <- rep_row_vector(1.0, 3) * (disc_mat2 .* block(wOthers[s], t-2, 1, 3, 4) .* reward_mat2 );
      }

      otherValue2[t+1,choice2[s,t]]   <- sum( othW2      .* otherReward_vec2 );
      otherValue2[t+1,3-choice2[s,t]] <- sum( (1-othW2 ) .* otherReward_vec2 );
      otherValue2[t+1] <- (1 ./ ( 1 + exp(-otherValue2[t+1]))) * 2 - 1;
      
    }  // trial loop
  }    // subject loop
}
