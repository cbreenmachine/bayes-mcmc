# time_wilcoxon_ranksum_test.R
# Can we sieve
library(data.table)
library(tidyverse)

# Number of CpGs we want
G <- 10

# Methylated and coverage
M <- fread("../dataRaw/chr17.M.bed") %>% arrange(chromStart)
Cov <- fread("../dataRaw/chr17.Cov.bed") %>% arrange(chromStart)

# Samplesheet
df <- read_csv("../dataRaw/masterSamplesheet.csv", show_col_types = F) %>%
  dplyr::filter(sample_id %in% colnames(M)) %>% distinct()

# Only use positions with decent coverage
medcov <- apply(Cov[ , -"chromStart"], MARGIN=1, FUN=median)
deep_index <- which(medcov >= 5)

# Load pvals and make sure the positions are lined up correctly
pvals <- fread("../dataRaw/pvals.bed") %>%
  dplyr::filter(chr == "chr17") %>%
  dplyr::mutate(chromStart = as.numeric(start)) %>%
  arrange(chromStart)

gc()

# Check order
all(pvals$start == M$chromStart)

# Now we can use positional indexing
sig_index <- which(pvals$lfdr < 0.05)

# Close neighbors measured as average
# neighbor distance
pos_mean_dist <- zoo::rollmean(diff(c(0, pvals$chromStart)), G)
close_index <- which(pos_mean_dist < 100)

# Loci where desqucning is deep and there are significant sites
best_index <- Reduce(intersect, list(deep_index, sig_index, close_index))
head(pvals[best_index, ])

# Pick one to center on
# 5100 is a reasonable choice...
# 5300 also...
anchor_ix <- best_index[5300] # arbitrarily pick one that satisfies conditions
keepix <- (anchor_ix-5):(anchor_ix+5) # indices of what to keep

# Probability matrix and get the columns in the right order
y <- t(M[keepix, -"chromStart"])
c <- t(Cov[keepix, -"chromStart"])

P.mat <- M[keepix, -"chromStart"] / Cov[keepix, -"chromStart"]
P.mat$position <- M$chromStart[keepix]

# Other data to throw into df
tmp.pvals <- dplyr::select(pvals, c("chromStart", "lfdr"))

df$bmi_std <- (df$bmi - mean(df$bmi, na.rm=T)) / sd(df$bmi, na.rm = T)
df$bmi_std[is.na(df$bmi_std)] <- 0

tmp.df <- dplyr::select(df, c("sample_id", "age_at_visit", "diagnostic_group", "sex", "bmi_std")) %>%
  transmute(sample_id = as.character(sample_id),
            is_load = ifelse(diagnostic_group == "LOAD", 1, 0),
            age_std = (age_at_visit - mean(age_at_visit)) / sd(age_at_visit),
            bmi_std,
            is_male = ifelse(sex == "Male", 1, 0))

# Data munging
data <- P.mat %>%
  pivot_longer(-position, names_to = "sample_id", values_to = "me.perc") %>%
  left_join(tmp.pvals, by = c("position"="chromStart")) %>%
  left_join(tmp.df, by = "sample_id")

p <- data %>%
  dplyr::mutate(Group = ifelse(is_load == 1, "Case", "Control")) %>%
  ggplot(aes(x = as.factor(position),
             y = me.perc,
             fill = Group)) +
  geom_boxplot(alpha = 0.9, outlier.size = 0.5) +
  theme_minimal() +
  ylab("Methylation (%)") +
  xlab("Genomic position (nt)") +
  ggsci::scale_fill_nejm() +
  ylim(c(0, 1)) +
  theme(axis.text.x = element_text(angle = 15, vjust = 0.5, hjust=0.35),
        plot.background = element_rect(color = "white"))

p
cowplot::save_plot(filename = "../figs/elevenLociBoxPlots.png", plot = p,
                   base_width = 8)

# Pull out the design matrix
# Data munging
X <- P.mat[1, ] %>%
  pivot_longer(-position, names_to = "sample_id", values_to = "me.perc") %>%
  left_join(tmp.pvals, by = c("position"="chromStart")) %>%
  left_join(tmp.df, by = "sample_id") %>%
  dplyr::select(is_load, age_std, bmi_std, is_male)

N <-nrow(X)

save(list = c("X", "G", "N", "c", "y"),
     file="../dataDerived/20221208-tenLoci-v2.RData")
