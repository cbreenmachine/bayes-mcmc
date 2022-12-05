library(rstan)
library(tidyverse)

load("../dataDerived/singleLocus.RData")

single_locus_model <-
  rstan::stan_model(file = "model-v1-oneCpG.stan")

# Assign variables (should be in its own script)
N <- nrow(single.locus.df)
K <- ncol(single.locus.df) - 3

X <- single.locus.df[ , c("isMale", "hasLOAD", "age_std")]
y <- single.locus.df$methylation

scale_alpha <- 1
scale_beta <- rep(1, K)
loc_sigma <- 2

model_args <- list(N = N, K = K,
                   X = X, y = y,
                   scale_alpha = scale_alpha,
                   scale_beta = scale_beta,
                   loc_sigma = loc_sigma)

# Let's say we want a total of 5000 samples at the end of the day
# We can run 4 chains for 2500 = 5000/2 iterations each and discard the
# first half of each chain as warmup.
N_samples <- 5000
fit <-
  rstan::sampling(object = single_locus_model,
                  data = model_args,
                  chains = 4, iter = N_samples/2)
