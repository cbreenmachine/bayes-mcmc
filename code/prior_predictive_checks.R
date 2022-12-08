library(tidyverse)
library(ggsci)

load("../dataDerived/20221130-tenLoci-v1.RData")

# Data munging for plotting
mm <- t(y / C)
colnames(mm) <-  1:ncol(mm)
mm <- mm %>%
  cbind(X) %>%
  as.data.frame() %>%
  pivot_longer(cols = -c(load_status, age_std, is_male),
               values_to = "m.perc",
               names_to = "locus") %>%
  dplyr::mutate(
    locus = factor(locus, levels = as.character(1:ncol(y))),
    Status = ifelse(load_status == 1, "Alz", "Ctrl"),
    Sex = ifelse(is_male == 1, "M", "F")
    )

# Plotting
p <- mm %>%
  group_by(locus) %>%
  ggplot(aes(x = locus, y = m.perc, fill = Status, color = Status)) +
  geom_boxplot() +
  theme_minimal() +
  scale_color_nejm() +
  scale_fill_nejm() +
  ylab("Methylation (%)") +
  xlab("Locus") +
  theme(legend.position = c(0.9, 0.15),
        plot.background = element_rect(fill = "white", color = "white")) +
  ggtitle("Ten CpG loci used in model development")

p

cowplot::save_plot(filename = "../figs/tenLociBoxPlots.png", p, base_width = 7)

# Need for age as covariate?
mm %>%
  ggplot(aes(x = age_std, y = m.perc)) +
  geom_point() +
  facet_wrap(.~locus) +
  theme_minimal()


# Plotting
p <- mm %>%
  ggplot(aes(x = locus, y = m.perc, fill = Sex, color = Sex)) +
  geom_boxplot() +
  theme_minimal() +
  ylab("Methylation (%)") +
  xlab("Locus") +
  theme(legend.position = c(0.9, 0.15),
        plot.background = element_rect(fill = "white", color = "white")) +
  ggtitle("Ten CpG loci used in model development")

mm %>%
  ggplot(aes(x = m.perc, fill = Sex)) +
  geom_histogram(position = "identity",
                 alpha = 0.5, bins = 50,
                 aes(y = after_stat(density))) +
  theme_minimal() +
  facet_wrap(.~locus)

mm %>%
  ggplot(aes(x = m.perc)) +
  geom_histogram(bins = 40, aes(y = after_stat(density)))

