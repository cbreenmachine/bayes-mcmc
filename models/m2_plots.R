library(tidyverse)
library(rstan)

load("m2_fit.RData")
# load("../dataDerived/20221208-tenLoci-v2.RData")
bayes.fit

check_hmc_diagnostics(bayes.fit)

phi_samples <-
  rstan::extract(bayes.fit, pars = paste0("phi[", 1:G, "]"), permuted = T) %>%
  as.data.frame()

colMeans(phi_samples)


# Spatial dependence (alpha) ----------------------------------------------
# extract data
alpha_samples <- extract(bayes.fit, pars = "alpha")[[1]]

# Plot
png("../figs/m2_alpha_hist.png", width = 6, height = 4, units = "in", res = 450)
hist(alpha_samples, prob = T,
     breaks = 50, main = "",
     xlab = expression(alpha))
abline(v = median(alpha_samples), col = "red", lwd = 2)
dev.off()


# Spatial smoothing parameters (phi)_j -------------------------------------

phis <- rstan::summary(bayes.fit, pars = paste0("phi[", 1:G, "]"))$summary

data.frame(phis) %>%
  mutate(j = 1:G) %>%
  ggplot(aes(x = j, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = `X2.5.`, ymax = `X97.5.`),
              fill = "red", alpha = 0.3) +
  theme_minimal() +
  xlab("")



# Beta LOAD ---------------------------------------------------------------

extract_beta <- function(ix, nn){
  z <- rstan::summary(bayes.fit, pars = paste0("beta[", ix, ",", 1:G, "]"))$summary

  # Select the mean, lower and upper bounds
  out <- data.frame(z) %>% dplyr::select(c(1, 4, 8))

  # Rename from the "cleaned" versions
  names(out) <- c("mean", "lower.95", "upper.95")
  out$locus <- 1:nrow(out)
  out$predictor <- nn
  out
}

betas.df <- rbind(
  extract_beta(1, "is_load"),
  extract_beta(2, "age_std"),
  extract_beta(3, "bmi_std"),
  extract_beta(4, "is_male")
)

pred.labs <- c("Age", "BMI", "AD status", "Sex")
names(pred.labs) <- c("age_std", "bmi_std", "is_load", "is_male")

p <- betas.df %>%
  ggplot(aes(x = locus, y = mean)) +
  geom_point() +
  geom_ribbon(aes(ymin = lower.95, ymax = upper.95),
              fill = "red", alpha = 0.3) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  ylab("Posterior mean") +
  xlab("Locus") +
  facet_wrap(.~predictor, scales = "free",
             labeller = labeller(predictor = pred.labs))


cowplot::save_plot(filename = "../figs/m2_predictor_bounds.png", p)

# Beta associated with LOAD -----------------------------------------------

me.probs <- y / c

fit_summary <- summary(bayes.fit)

me.probs.hat <-
  as.data.frame(fit_summary) %>%
  rownames_to_column("parameter") %>%
  dplyr::filter(str_detect(parameter, "mu")) %>%
  dplyr::mutate(parameter = str_remove(parameter, "mu\\[")) %>%
  dplyr::mutate(parameter = str_remove(parameter, "\\]")) %>%
  separate(parameter, ",", into = c("sample", "locus")) %>%
  dplyr::select(c("sample", "locus", "summary.mean")) %>%
  pivot_wider(id_cols = sample, names_from = locus, values_from = summary.mean) %>%
  select(-sample) %>%
  as.matrix() %>%
  boot::inv.logit()


xx <- as.vector(me.probs)
yy <- as.vector(me.probs.hat)

rho <- round(cor(xx, yy), 2)

plot(xx, yy, pch = 20,
     xlab = "Ground truth methylation",
     ylab = "Predicted methylation",
     xlim = c(0,1), ylim = c(0, 1),
     main = "Title")
