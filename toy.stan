data {
    int N;
    vector[N] x;
    vector[N] y;  
}

parameters {
    real beta_mu;
    real<lower=0> beta_sd;
    vector[11] beta_raw;
}

transformed parameters {
    vector[11] beta;
    beta <- beta_mu + beta_sd * beta_raw;
}

model {
    beta_mu  ~ normal(0,1);
    beta_sd  ~ cauchy(0,5);
    beta_raw ~ normal(0,1);
    
    y ~ normal(beta[1]+beta[2]*x+beta[3]*log(x), beta[4]);
}