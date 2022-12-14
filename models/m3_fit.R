library(rstan)

N_samples <- 10000
N_cores <- 5

ifile <- "../dataDerived/20221208-tenLoci-v2.RData"
ofile <- "m2_fit.RData"

# RUN THIS BLOCK AT ONCE
load(ifile)
P <- ncol(X)

G <- 11

W <- matrix(data = 0, nrow = G, ncol = G)

for (i in 1:(G-1)){
  W[i, i+1] <- 1
  W[i+1, i] <- 1
}

W

# Set parameters and data
args <- list(N = N, G = G, P = P,
             W = W,
             X = X, y = y, c = c)

# SAMPLE

model <- rstan::stan_model("m2_binomial_covariates_re.stan")

start_time <- Sys.time()
bayes.fit <-
  rstan::sampling(object = model,
                  data = args,
                  chains = N_cores,
                  iter = N_samples,
                  cores = N_cores)
elapsed_time <- Sys.time() - start_time


save(list = ls(), file = ofile)
