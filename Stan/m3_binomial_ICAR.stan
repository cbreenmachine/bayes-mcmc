//
// Capture the over-dispersion in the means

// LOAD random effect
// CpG specific random effect

data {
  int<lower=1> N; // number of subjects
  int<lower=1> G; // number of CpGs
  int<lower=1> P; // number of predictors (don't include intercept)

  matrix[N, P] X; // covariate --> LOAD (1) and control (0)

  matrix<lower=0, upper=1>[G, G] W;

  int y[N, G]; // count outcomes
  int c[N, G]; // (known) number of trials in Binomial
}

transformed data {
  // used in the MVN prior later on
  vector[G] zeros;

  // Count up the number of neighbors
  // each CpG locus has. It will be 1 or 2
  // until something changes
  matrix<lower=0>[G, G] D;
  {
    vector[G] W_rowsums; //
    for (i in 1:G) {
      W_rowsums[i] = sum(W[i, ]);
    }
    D = diag_matrix(W_rowsums);
  }
  zeros = rep_vector(0, G);
}


parameters {
  matrix[P, G] beta; // one length P vector of betas for each CpG
  vector[G] phi; //spatial component for each CpG
  // real<lower=0.9, upper=1>alpha; // uniform prior on dependence param

  real<lower=0> tau; // precision of spatial RV
}


transformed parameters {
  matrix[N, G] mu; // adjusted mean

  for (i in 1:N){
    for (j in 1:G){
      mu[i,j] = X[i, ]*beta[ ,j] + phi[j];
    }
  }
}


model {
  // For some reason, setting alpha==1 does not work
  phi ~ multi_normal_prec(zeros, tau * (D - 0.99*W));

  // cumbersome, can we vectorize?
  for (p in 1:P){
    for (j in 1:G){
      beta[p,j] ~ normal(0, 1);
    }
  }

  // justification
  tau ~ gamma(2, 2);

  for (i in 1:N){
    for (j in 1:G){
      y[i,j] ~ binomial_logit(c[i,j], mu[i,j]);
    }
  }
}
