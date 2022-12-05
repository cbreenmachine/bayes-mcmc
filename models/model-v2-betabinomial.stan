//Defines a very simple one CpG model
//Use

data {
  int<lower=0> N; // Number of participants
  int m[N]; // methylation percentage

  // Design-related
  int<lower=0> K; // number of covariates (not including intercept)
  matrix [N, K] X; // design matirx (no intercept)
  int C[N];

}

parameters {
  vector[K] b; // including the intercept
  real phi; // dispersion
  real b0; // intercept
}

transformed parameters {

  vector[N] mu;
  vector[N] alpha; //shape 1
  vector[N] beta; //shape 2

  for (n in 1:N)
    mu[n] = inv_logit(b0 + X[n, ]*b); //logit link

  alpha = mu * phi;
  beta = (1-mu) * phi;
}


model {
  m ~ beta_binomial(C, alpha, beta);
}

generated quantities {
  vector[N] m_rep;
  for (n in 1:N){
    m_rep[n] = beta_binomial_rng(C[n], alpha[n], beta[n]);
  }
}

