library(tidyverse)
library(latex2exp)


samples.df <- read_csv("../dataRaw/masterSamplesheet.csv", show_col_types = F)
covariates.df <- samples.df %>%
  dplyr::mutate(isMale = ifelse(sex == "Male", 1, 0),
                hasLOAD = ifelse(diagnostic_group == "LOAD", 1, 0),
                age = age_at_visit,
                sample_id = as.character(sample_id)) %>%
  dplyr::select(c(sample_id, age, isMale, hasLOAD))

head(covariates.df)


df <- read_csv("../dataDerived/modelInputData.csv", show_col_types = F) %>%
  pivot_longer(cols = -position, values_to="methylation", names_to="sample_id") %>%
  left_join(covariates.df, by = "sample_id")


var.df <- df %>%
  group_by(position) %>%
  summarize(V = var(methylation))

vv <- var.df$V

vv.rank <-rank(-vv)
vv[vv.rank == 12]

pos <- var.df$position[vv.rank == 12]


single.locus.df <- df %>%
  dplyr::filter(position == pos) %>%
  dplyr::transmute(position, sample_id, methylation, isMale, hasLOAD,
                   age_std = (age - mean(age)) / sd(age))


save(single.locus.df, file="../dataDerived/singleLocus.RData")

png("../figs/distVariances.png")
hist(vv, breaks=50,
     probability = T,
     main = "Estimated variance in methylation by locus",
     xlab = "Estimated variance")
dev.off()


summary.df <- df %>%
  group_by(position, hasLOAD) %>%
  summarize(mean.me = mean(methylation),
            sd.me = sd(methylation)) %>%
  mutate(low = mean.me - 1.96 * sd.me,
         high = mean.me + 1.96 * sd.me)




summary.df[900:1000, ] %>%
  ggplot(aes(x = position, y = mean.me, ymin = low, ymax = high, color = as.factor(hasLOAD))) +
  geom_point() +
  geom_errorbar(alpha = 0.7) +
  theme_minimal() +
  xlab("Position (index)") +
  ylab(TeX("Mean methylation $\\pm 1.96 \\times SD$"))
