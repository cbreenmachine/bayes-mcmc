library(tidyverse)
library(rstan)

# ifile <- "../ModelFits/m1_binomial_result.RData"
# ofile_prefix <- "../Figures/m1_"

ifile <- "../ModelFits/m2_binomial_CAR_result.RData"
ofile_prefix <- "../Figures/m2_"

# ifile <- "../ModelFits/m3_binomial_ICAR_result.RData"
# ofile_prefix <- "../Figures/m3_"

print(ifile)
load(ifile)

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
  facet_wrap(.~predictor, scales = "free", nrow = 1,
             labeller = labeller(predictor = pred.labs))


cowplot::save_plot(filename = paste0(ofile_prefix, "beta_with_ci.png"), p, base_width = 10)

# Raw vs adjusted methylation % -----------------------------------------------

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

p <- data.frame(raw = xx, y.hat = yy) %>%
  ggplot(aes(x = raw, y = y.hat)) +
  geom_point() +
  theme_minimal() +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  xlab("Ground truth methylation (%)") +
  ylab("Predicted methylation (%)")

cowplot::save_plot(filename = paste0(ofile_prefix, "ground_truth_vs_adjusted.png"), p)


# Trace plot for beta[6] i.e. LOAD ----------------------------------------
# rstan::traceplot(bayes.fit, pars = "phi[6]")

png(paste0(ofile_prefix, "rhat.png"))
stan_rhat(bayes.fit)
dev.off()

png(paste0(ofile_prefix, "ess.png"))
stan_ess(bayes.fit)
dev.off()
