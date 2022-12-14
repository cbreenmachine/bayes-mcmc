library(rstan)

N_samples <- 10000
N_cores <- 5

ifile <- "../dataDerived/20221208-tenLoci-v2.RData"
ofile <- "m1_fit.RData"

# RUN THIS BLOCK AT ONCE
load(ifile)
P <- ncol(X)
G <- 11

# Set parameters and data
args <- list(N = N, G= G, P =P,
             X = X, y = y, c = c)

# SAMPLE
model <- rstan::stan_model("m1_binomial.stan")

start_time <- Sys.time()
bayes.fit <-
  rstan::sampling(object = model,
                  data = args,
                  chains = N_cores,
                  iter = N_samples,
                  cores = N_cores)
elapsed_time <- Sys.time() - start_time

save(list = ls(), file = ofile)
