library(rstan)

N_samples <- 5000
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

bayes.fit <-
  rstan::sampling(object = model,
                  data = args,
                  chains = , iter = N_samples/2,
                  cores = 6)

bayes.fit

pairs(bayes.fit, pars = c("beta[1,1]", "beta[1,2]"))

load_effects <- as.data.frame(bayes.fit) %>%
  dplyr::select(starts_with("beta[1,"))

quantile(load_effects[ ,2], c(0.05, 0.95))

hist(load_effects[,2])

save(list = ls(), file = ofile)
