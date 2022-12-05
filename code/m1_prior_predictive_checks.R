
# Load and reshape
load("../dataDerived/20221130-tenLoci-v1.RData")
C <- t(C); y <- t(y)

# Constants
N <- nrow(y) # number of subjects
G <- ncol(y) # number of loci

# Methylation probability
p <- unlist(y / C)
hist(p)

# Priors...
mu_sample <- runif(n = N*G, min = 0, max=1)
# mu_sample <- rnorm(n = N*G, mean = 0.75, sd = 0.2)
lambda_sample <- EnvStats::rpareto(n = N*G, 1, 1.5)
# lambda_sample <- runif(n = N*G, min=0.01, max=0.5)
# lambda_sample <- 20

alpha <- mu_sample * lambda_sample
beta <- (1 - mu_sample) * lambda_sample

p_sample <- rbeta(n = length(alpha), shape1 = alpha, shape2 = beta)
y_sample <- rbinom(n = length(p_sample), prob = p_sample, size = unlist(C))

hist(y_sample / C)
