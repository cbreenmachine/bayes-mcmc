library(rstan)
source("constants.R")
ofile <- "../ModelFits/m1_binomial_result.RData"

load(ifile)
P <- ncol(X)
G <- 11

# Set parameters and data
args <- list(N = N, G = G, P = P,
             X = X, y = y, c = c)

# SAMPLE
model <- rstan::stan_model("m1_binomial.stan")

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

#END
