// 10 CpG by 281 particiapnt dataset

data {
  int<lower=0> N; // number of samples
  int<lower=0> G; // number of CpGs (likely fewer than 10)
  int<lower=0> P; // number of parameters in design matrix (baseline is 3)

  // int<lower=0> N_edges;
  // int<lower=0, upper=N> node1[N_edges]; // node1[i] adjacent to node2[i]
  // int<lower=0, upper=N> node2[N_edges]; // node1[i] adjacent to node2[i]

  matrix[N, P] X; // design matirx (no intercept)
  int C[N, G];
  int y[N, G];
}

parameters {
  real b0[G]; // intercepts for G CpGs
  matrix[P, G] b; // slopes for age, case-control, sex

  real<lower=0> psi; // dispersion common to all CpGs...
}

transformed parameters {
  // All CpG and patient-specific parameters
  matrix[N, G] alpha; //shape 1
  matrix[N, G] beta; //shape 2
  matrix<lower=0, upper=1>[N, G] mu; // mean

  for (i in 1:N) { // thru sample
    for (j in 1:G){ // thru CpG
      // print("loop iteration: ", i, " ", j);
      mu[i, j] = inv_logit(b0[j] + X[i, ]*b[ ,j]); //logit link
    }
  }

  // cast into shape1 and shape2
  alpha = mu * psi;
  beta = (1-mu) * psi;
}

model {

  for (i in 1:N) {
    for (j in 1:P){
      y[i, j] ~ beta_binomial(C[i,j], alpha[i,j], beta[i,j]);
    }
  }

}
