library(rstan)
library(tidyverse)

load("../dataDerived/singleLocus.RData")

single_locus_model <-
  rstan::stan_model(file = "model-v2-betabinomial.stan")

# Assign variables (should be in its own script)
N <- nrow(locus.df)
K <- ncol(locus.df) - 3 # no intercept included

X <- locus.df[ , c("hasLOAD", "isMale", "age_std")]
m <- locus.df$M
C <- locus.df$C


model_args <- list(N = N, K = K,
                   X = X, m = m,
                   C = C)

# Let's say we want a total of 5000 samples at the end of the day
# We can run 4 chains for 2500 = 5000/2 iterations each and discard the
# first half of each chain as warmup.
N_samples <- 25000
bayes.fit <-
  rstan::sampling(object = single_locus_model,
                  data = model_args,
                  chains = 4, iter = N_samples/2,
                  cores = 4)


pairs(fit, pars = c("b0", "b[1]", "b[2]", "b[3]", "phi"))

hist(as.data.frame(fit)$`m_rep[100]`)
