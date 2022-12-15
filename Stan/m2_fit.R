library(rstan)

source("constants.R")
ofile <- "../ModelFits/m2_binomial_CAR_result.RData"

# RUN THIS BLOCK AT ONCE
load(ifile)
P <- ncol(X)
G <- 11

# CReate weight matrix
W <- matrix(data = 0, nrow = G, ncol = G)

for (i in 1:(G-1)){
  W[i, i+1] <- 1
  W[i+1, i] <- 1
}

# Set parameters and data
args <- list(N = N, G = G, P = P,
             W = W,
             X = X, y = y, c = c)

# SAMPLE
model <- rstan::stan_model("m2_binomial_CAR.stan")

start_time <- Sys.time()
bayes.fit <-
  rstan::sampling(object = model,
                  data = args,
                  chains = N_cores,
                  iter = N_samples,
                  cores = N_cores,
                  seed = seed)
elapsed_time <- Sys.time() - start_time
save(list = ls(), file = ofile)

# END
