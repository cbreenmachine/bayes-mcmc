library(tidyverse)

# Generate N row by G column
G  <- 25 # number of CpG loci
N <- 100 # number of samples (N/2 cases and controls)

# Realized coverage (nmber of trials)
c <- matrix(data = 30, nrow = N, ncol = G)

# Average methylation per site
xx <- seq(from = -pi/2, to = pi/2, length.out = G)
baseline_pi <- cos(xx)



# (S1) DMR ----------------------------------------------------------------

# Create hypomethylated region
X <- matrix(rep(c(0,1), N/2), ncol=1)

p <- matrix(NA, nrow = N, ncol = G)
y <- matrix(NA, nrow = N, ncol = G)

set.seed(919)

for (i in 1:N){
  for (j in 1:G){
    z <- baseline_pi[j] + rnorm(1, mean=0, sd=0.05)

    # Only want the hypo methylation to hit the middle positions
    if (j > 9 & j < 16) {
      z2 <- z + X[i]*rnorm(n=1, mean= -0.2, sd=0.05)
    } else {
      z2 <- z - 0.4
    }
    p[i,j] <- z2
  }
}

p[p < 0] <- 0.01
p[p > 1] <- 0.99

for (i in 1:N){
  for (j in 1:G){
    y[i,j] <- rbinom(n = 1, size = c[i,j], prob = p[i,j])
  }
}

# OUTPUT
save(list = c("N", "G", "X", "y", "c"),
     file = "../DataSimulated/s1.RData")


# Plot to make sure this looks sensible
dat <- data.frame(cbind(p, X))
colnames(dat) <- c(1:G, "group")

# Data plot
p1 <- dat %>%
  pivot_longer(-group) %>%
  mutate(x = as.numeric(name),
         Group = as.factor(group)) %>%
  ggplot(aes(x = x, y = value, color = Group)) +
  geom_point() +
  theme_minimal() +
  xlab("Locus") +
  ylab("Simulated p") +
  ggsci::scale_color_nejm() +
  ggtitle("(S1) Diffentially methylated region") +
  theme(legend.position = c(0.8, 0.9))


p1



# (S2) Switch ----------------------------------------------------------------

# Shift the region down a bit
baseline_pi <- 0.5*cos(xx)
# baseline_pi[10:15] <- baseline_pi[10:15] - 0.5

p <- matrix(NA, nrow = N, ncol = G)
y <- matrix(NA, nrow = N, ncol = G)

set.seed(919)
for (i in 1:N){
  for (j in 1:G){
    z <- baseline_pi[j] + rnorm(1, mean=0, sd=0.05)

    # Only want the hypo methylation to hit the middle positions
    if (j > 8 & j < 12) {
      z2 <- z - X[i]*rnorm(n=1, mean = 0.2, sd=0.05)
    } else if (j > 12 & j < 16){
      z2 <- z + X[i]*rnorm(n=1, mean= 0.2, sd=0.05)
    } else{
      z2 <- z
    }

    p[i,j] <- z2
  }
}

p[p < 0] <- 0.01
p[p > 1] <- 0.99



for (i in 1:N){
  for (j in 1:G){
    y[i,j] <- rbinom(n = 1, size = c[i,j], prob = p[i,j])
  }
}


# OUTPUT
save(list = c("N", "G", "X", "y", "c"),
     file = "../DataSimulated/s2.RData")



# Plot to make sure this looks sensible
dat <- data.frame(cbind(p, X))
colnames(dat) <- c(1:G, "group")

# Data plot
p2 <- dat %>%
  pivot_longer(-group) %>%
  mutate(x = as.numeric(name),
         Group = as.factor(group)) %>%
  ggplot(aes(x = x, y = value, color = Group)) +
  geom_point() +
  theme_minimal() +
  xlab("Locus") +
  ylab("Simulated p") +
  ggsci::scale_color_nejm() +
  ggtitle("(S2) Switching event") +
  theme(legend.position = c(0.8, 0.9))

p2


# Combine and output ------------------------------------------------------

z <- cowplot::plot_grid(p1, p2)
cowplot::save_plot(filename = "../Figures/simulated_points.png", z, base_width = 8)
