//
// Capture the over-dispersion in the means

// LOAD random effect
// CpG specific random effect

data {
  int<lower=0> N; // number of subjects
  int<lower=0> G; // number of CpGs
  int<lower=0> P; // number of predictors (don't include intercept)

  matrix[N, P] X; // covariate --> LOAD (1) and control (0)

  int y[N, G]; // count outcomes
  int c[N, G]; // (known) number of trials in Binomial
  real<lower=0> tau_std; // standard deviation of random effect with mean 0
}

parameters {
  vector[G] alpha; // intercept
  matrix[P, G] beta; // slope (LOAD effect)
  real tau; // random effect for each CpG

}

transformed parameters {
  matrix[N, G] theta; // adjusted mean

  for (i in 1:N){
    for (j in 1:G){
      theta[i,j] = alpha[j] + X[i, ]*beta[ ,j] + tau;
    }
  }
}

model {
  // prior specification
  tau ~ normal(0, tau_std);

  for (i in 1:N){
    for (j in 1:G){
      y[i,j] ~ binomial_logit(c[i,j], theta[i,j]);
    }
  }
}
