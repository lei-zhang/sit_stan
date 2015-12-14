data {
    int N;
    vector[N] x;
    vector[N] y;
}

transformed data {
    int B;

    B <- 2;
}

parameters {
    real beta_mu;
    real<lower=0> beta_sd;
    vector[B] beta_raw;
    matrix[100,4] A[129];
}

transformed parameters {
    vector[2] beta;
    beta <- beta_mu + beta_sd * beta_raw;
}

model {
    beta_mu  ~ normal(0,1);
    beta_sd  ~ cauchy(0,5);
    beta_raw ~ normal(0,1);
    
    y ~ normal(beta[1]+beta[2]*x, 2);
}