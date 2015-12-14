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
  matrix<lower=0,upper=1>[nTrials,4]  otherWith2[nSubjects];
  matrix<lower=0.25,upper=1>[nTrials,4]  wOthers[nSubjects]; 
}

transformed data {
  vector[2] initV;    // initial values for V
  int B;              // number of predictors
  vector[3] pwr;      // power
  
  initV  <- rep_vector(0.0,2); 
  B      <- 6;
  pwr[1] <- 2;
  pwr[2] <- 1;
  pwr[3] <- 0;
}

parameters {
  // subject-level raw parameters, follows norm(0,1), for later Matt Trick
  vector<lower=0,upper=1>[nSubjects] lr;
  vector<lower=0,upper=1+1e-10>[nSubjects] disc;
  vector[nSubjects] beta[B];
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
  row_vector[4] othW;

  for (i in 1:B) {
    beta[i] ~ normal(0,5);
  } 
  
  
  // subject loop and trial loop
  for (s in 1:nSubjects) {
    myValue[1]    <- initV;
    otherValue[1] <- initV;

    for (t in 1:nTrials) {
      valfun1 <- beta[1,s]*myValue[t] + beta[2,s]*otherValue[t];
      choice1[s,t] ~ categorical_logit( valfun1 );

      valdiff <- myValue[t,choice1[s,t]] - myValue[t,3-choice1[s,t]];
      valfun2 <- beta[3,s] + beta[4,s]*valdiff + beta[5,s]*wgtWith[s,t] + beta[6,s]*wgtAgst[s,t];
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
        otherValue[t+1,choice2[s,t]]   <- sum( wOthers[s,t] .* othW .* otherReward2[s,t] );
        otherValue[t+1,3-choice2[s,t]] <- sum( wOthers[s,t] .* (1-othW) .* otherReward2[s,t] );
      } else if (t==2) {
        otherValue[t+1,choice2[s,t]]   <- sum( wOthers[s,t] .* othW .* otherReward2[s,t]     + wOthers[s,t-1] .* othW .* otherReward2[s,t-1]*disc[s] );
        otherValue[t+1,3-choice2[s,t]] <- sum( wOthers[s,t] .* (1-othW) .* otherReward2[s,t] + wOthers[s,t-1] .* (1-othW) .* otherReward2[s,t-1]*disc[s] );
      } else {
        disc_mat <- rep_matrix(exp(log(disc[s])*pwr),4); // replicate by 4 columns 
        otherValue[t+1,choice2[s,t]]   <- sum( disc_mat .* block(wOthers[s], t-2, 1, 3, 4) .* rep_matrix(othW,3)   .* block(otherReward2[s], t-2, 1, 3, 4) );
        otherValue[t+1,3-choice2[s,t]] <- sum( disc_mat .* block(wOthers[s], t-2, 1, 3, 4) .* rep_matrix(1-othW,3) .* block(otherReward2[s], t-2, 1, 3, 4) );
      }
    }  // trial loop
  }    // subject loop
}

generated quantities {
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
  row_vector[4] othW2;
  int<lower=0,upper=1> c_rep[nSubjects, nTrials];

  for (s in 1:nSubjects) {
    myValue2[1]    <- initV;
    otherValue2[1] <- initV;
    log_likc1[s]   <- 0;
    log_likc2[s]   <- 0;

    for (t in 1:nTrials) {
      valfun1_gen  <- beta[1,s]*myValue2[t] + beta[2,s]*otherValue2[t];
      log_likc1[s] <- log_likc1[s] + categorical_logit_log(choice1[s,t], valfun1_gen);

      valdiff_gen  <- myValue2[t,choice1[s,t]] - myValue2[t,3-choice1[s,t]];
     
      valfun2_gen  <- beta[3,s] + beta[4,s]*valdiff_gen + beta[5,s]*wgtWith[s,t] + beta[6,s]*wgtAgst[s,t];
      log_likc2[s] <- log_likc2[s] + bernoulli_logit_log(chswtch[s,t], valfun2_gen);
      
      c_rep[s,t]   <- bernoulli_rng( inv_logit(valfun2_gen) );

      pe2[t]   <-  reward[s,t] - myValue2[t,choice2[s,t]];
      penc2[t] <- -reward[s,t] - myValue2[t,3-choice2[s,t]];
      
      myValue2[t+1,choice2[s,t]]   <- myValue2[t,choice2[s,t]]   + lr[s] * pe2[t];
      myValue2[t+1,3-choice2[s,t]] <- myValue2[t,3-choice2[s,t]] + lr[s] * penc2[t];

      othW2 <- otherWith2[s,t]; 

      if (t==1) {
        otherValue2[t+1,choice2[s,t]]   <- sum( wOthers[s,t] .* othW2 .* otherReward2[s,t] );
        otherValue2[t+1,3-choice2[s,t]] <- sum( wOthers[s,t] .* (1-othW2) .* otherReward2[s,t] );
      } else if (t==2) {
        otherValue2[t+1,choice2[s,t]]   <- sum( wOthers[s,t] .* othW2 .* otherReward2[s,t]     + wOthers[s,t-1] .* othW2 .* otherReward2[s,t-1]*disc[s] );
        otherValue2[t+1,3-choice2[s,t]] <- sum( wOthers[s,t] .* (1-othW2) .* otherReward2[s,t] + wOthers[s,t-1] .* (1-othW2) .* otherReward2[s,t-1]*disc[s] );
      } else {
        disc_mat2 <- rep_matrix(exp(log(disc[s])*pwr),4);
        otherValue2[t+1,choice2[s,t]]   <- sum( disc_mat2 .* block(wOthers[s], t-2, 1, 3, 4) .* rep_matrix(othW2,3)   .* block(otherReward2[s], t-2, 1, 3, 4) );
        otherValue2[t+1,3-choice2[s,t]] <- sum( disc_mat2 .* block(wOthers[s], t-2, 1, 3, 4) .* rep_matrix(1-othW2,3) .* block(otherReward2[s], t-2, 1, 3, 4) );
      }
    }  // trial loop
  }    // subject loop
}