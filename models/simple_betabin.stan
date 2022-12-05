// LOAD random effect
// CpG specific random effect

data {
  int<lower=0> N; // number of subjects
  int<lower=0> G; // number of CpGs
  int X[N]; // covariate --> LOAD (1) and control (0)

  int y[N, G]; // count outcomes
  int C[N, G]; // (known) number of trials in Binomial

  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]
}

parameters {
  vector[G] beta0; // intercept
  vector[G] beta1; // slope

  real<lower=0> nu; // dispersion, for now keep constant
  real<lower=0> tau_theta;   // precision of heterogeneous effects
  real<lower=0> tau_phi;     // precision of spatial effects

  vector[G] phi;         // spatial effects
  vector[G] theta;  // heterogeneous effects

}

transformed parameters {

  real<lower=0> sigma_theta = inv(sqrt(tau_theta));  // convert precision to sigma
  real<lower=0> sigma_phi = inv(sqrt(tau_phi));      // convert precision to sigma

  matrix[N, G] mu; // mean
  matrix[N, G] alpha; // shape 1
  matrix[N, G] beta; // shape 2

  for (n in 1:N){
    for (g in 1:G){
      mu[n, g] = inv_logit(beta0[g] + beta1[g]*X[n] + phi[g]*sigma_phi + theta[g]*sigma_theta);
      alpha[n, g] = mu[n, g] * nu;
      beta[n, g] = (1-mu[n, g]) * nu;
    }
  }
}

model {

  theta ~ normal(0, 5);

  for (n in 1:N) {
    for (g in 1:G){
      y[n, g] ~ beta_binomial(C[n, g], alpha[n, g], beta[n, g]);
    }
  }

  // NOTE:  no prior on phi_raw, it is used to construct phi
  // the following computes the prior on phi on the unit scale with sd = 1
  target += -0.5 * dot_self(phi[node1] - phi[node2]);
  // soft sum-to-zero constraint on phi)
  sum(phi) ~ normal(0, 0.001 * N);  // equivalent to mean(phi) ~ normal(0,0.001)

  theta ~ normal(0, 1);
  tau_theta ~ gamma(3.2761, 1.81);  // Carlin WinBUGS priors
  tau_phi ~ gamma(1, 1);            // Carlin WinBUGS priors
}
//
// generated quantities {
//   matrix[N, G] mu_q;
//
//   for (i in 1:N){ // thru sample
//     for (g in 1:G){ // thru CpG
//       mu_q[i,g] = inv_logit(b0[g] + b1[g]*x[i]);
//     }
//   }
// }
