library(rstan)

load("../dataDerived/20221130-tenLoci-v1.RData")

G <- nrow(C)
N <- ncol(C)

node1 <- 1:(G-1)
node2 <- 2:G

N_edges <- length(node1)

args <-
  list( C = t(C),
        y = t(y),
        X = X[ ,1],
        node1 = node1,
        node2 = node2,
        N_edges = N_edges,
        N = N, G = G)


# Let's say we want a total of 5000 samples at the end of the day
# We can run 4 chains for 2500 = 5000/2 iterations each and discard the
# first half of each chain as warmup.
N_samples <- 5000

model <- rstan::stan_model("simple_betabin.stan")

bayes.fit <-
  rstan::sampling(object = model,
                  data = args,
                  chains = 6, iter = N_samples/2,
                  cores = 6)

bayes.fit

pairs(bayes.fit, )
