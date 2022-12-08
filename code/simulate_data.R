# simultate_data.R
# Three situations
# 1. No differential effect
# 2.

library(tidyverse)

# Helper functions
logit <- function(z){log(z / (1-z))}
inv_logit <- function(z){exp(z) / (1 + exp(z))}

anchor <- 0.9
intended_diff <- c(-0.05, -0.1, -0.15, -0.2, -0.25,
                   -0.3, -0.35, -0.4, -0.45, -0.5)

C <- 30 # Universal coverage (# trials in Binomial experiment)
N <- 250 # Number of samples / subjects--half will be case,half control
G <- 10 # Number of CpGs

# Mean of Gaussian from which to sample for alpha (intercept)
# inv_logit(0) == 0.5 --> center at half methylated
alpha_norm_mean <- rep(0, G)

# Range of case-control effects. These move logit()
# away from the 0.5 anchor point
beta_norm_mean <- seq(-1.5, 1.5, length.out = G) # No effect

# Common to all simulations
X <- matrix(rep(NA, (G+1)*N), nrow = N)
colnames(X) <- c("x", paste0("locus", 1:G))
X[ ,1] <- rep(c(0,1), N/2) # case, control

# Sample data for each subject one by one
set.seed(919) # "919 MF" -Petey Pablo
for (n in 1:N){

  # Sample the intercept and slope (in this case, case control effect)
  alpha <- alpha_norm_mean + rnorm(G, mean=0, sd=0.1)
  beta <- beta_norm_mean # + rnorm(G, mean=0, sd = 0.1)

  # Generating process
  z <- alpha + beta * X[n, 1]
  theta <- inv_logit(z)
  X[n, 2:(G+1)] <- rbinom(n = G, size = rep(C, G), prob = theta)
}


data <- as.data.frame(X) %>%
  pivot_longer(cols = starts_with("locus"),
               values_to = "y",
               names_to = "Locus") %>%
  dplyr::mutate(x = as.factor(x),
                Locus = factor(Locus, levels = paste0("locus", 1:10)))


data %>%
  ggplot(aes(x = Locus, y = y / C, fill = x, color = x)) +
  geom_boxplot() +
  theme_minimal()

p
