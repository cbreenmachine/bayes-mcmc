library(tidyverse)
library()

load("../dataDerived/singleLocus.RData")

X <- cbind(1, locus.df[ ,c("hasLOAD", "isMale", "age_std")])
y <- locus.df$M / locus.df$C

data <- cbind(X, y)

fit <- glm('y ~ hasLOAD + isMale + age_std',
            data = data,
            family = quasibinomial(link = "logit"))


mu_psi_to_alpha_beta <- function(mu, phi){
  # helper to go from mean/dispersion to canoncial parameterization
  return(c(
    mu * phi,
    (1-mu) * phi
  ))
}


mu_psi_to_alpha_beta(mu, phi)





sample_from_prior <- function(mu, phi, X, b, C){

  # Parameters for beta
  params <- mu_psi_to_alpha_beta(mu, phi)
  alpha <- params[1]
  beta <- params[2]

}


alpha <- 9
beta <- 11

p.sample <- rbeta(n = nrow(locus.df), shape1 = alpha, shape2 = beta)
m.sample <- rbinom(n = 1000, size = locus.df$C, prob = p.sample)


# Plot the spread of methylation percentages at this site
hist(locus.df$M,
     breaks = 40,
     prob = T,
     xlab = "Methylation counts",
     main = "Methylation at locus A",
     col = rgb(1,1,0,1/4))

hist(m.sample, add=T, prob = T, breaks=40, col = rgb(1,0,0,1/4))

