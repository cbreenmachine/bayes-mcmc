library(rstan)

ifile <- "../dataDerived/20221130-tenLoci-v1.RData"
ofile <- "m2_fit.RData"

tau_std <- 0.1

# RUN THIS BLOCK AT ONCE
load(ifile)
y <- t(y); c <- t(C)
N <- nrow(y); G <- ncol(y)
X <- X[ ,1:3]
P <- ncol(X)

# Set parameters and data
args <- list(N = N, G= G, P =P,
             X = X, y = y, c = c,
             tau_std = tau_std)

# SAMPLE
N_samples <- 5000

model <- rstan::stan_model("m2_binomial_covariates_re.stan")

bayes.fit <-
  rstan::sampling(object = model,
                  data = args,
                  chains = 6, iter = N_samples/2,
                  cores = 6)

bayes.fit

pairs(bayes.fit, pars = c("beta[1,1]", "beta[1,2]",
                          "beta[1,3]", "beta[1,4]",
                          "beta[1,5]", "beta[1,6]"))



load_effects <- as.data.frame(bayes.fit) %>%
  dplyr::select(starts_with("beta[1,"))

quantile(load_effects[ ,2], c(0.05, 0.95))

hist(load_effects[,2])

save(list = ls(), file = ofile)
