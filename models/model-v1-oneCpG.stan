//Defines a very simple one CpG model
//Use

data {
  int<lower=0> N; // Number of participants
  vector[N] y; // methylation percentage

  // Design-related
  int<lower=0> K; // number of covariates (not including intercept)
  matrix [N, K] X; // design matirx (no intercept)

  // Prior variances
  real scale_alpha; // prior on alpha
  vector[K] scale_beta; // prior on spread of beta
  real loc_sigma;

}

parameters {
  real alpha;
  vector[K] beta;
  real sigma;
}

transformed parameters {
  vector[N] mu;
  mu = alpha + X * beta;
}


model {
  // priors
  alpha ~ normal(0., scale_alpha);
  beta ~ normal(0., scale_beta);

  // LIKELIHOOD
  y ~ bernoulli_logit(mu);
}
