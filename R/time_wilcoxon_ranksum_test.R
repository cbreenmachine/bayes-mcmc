# time_wilcoxon_ranksum_test.R
# Can we sieve
library(data.table)
library(tidyverse)

N <- 100000

M <- fread("../dataRaw/chr17.M.bed")
Cov <- fread("../dataRaw/chr17.Cov.bed")

df <- read_csv("../dataRaw/masterSamplesheet.csv", show_col_types = F) %>%
  dplyr::filter(sample_id %in% colnames(M)) %>% distinct()

# Probability matrix and get the columns in the right order
P <- M[1:N,-"chromStart"] / Cov[1:N, -"chromStart"]
cols <- as.character(df$sample_id)
P <- P[ , ..cols] # puts P in the same order as df

all(colnames(P) == df$sample_id)
P <- as.matrix(P)

load.ix <- df$diagnostic_group == "LOAD"
cont.ix <- df$diagnostic_group == "CONTROL"

wrapper <- function(i){
  xx <- P[i, load.ix]; yy <- P[i, cont.ix]
  wilcox.test(xx, yy)$p.val
}

wrs.p <- do.call(rbind, parallel::mclapply(X = 1:N, FUN = wrapper, mc.cores = 6))

my_lab <- paste("Distribution of Wilcoxon\nRank-Sum pvals at", nrow(P), "loci")

tiff("../figs/wilcoxonPvals.tiff", width = 6, height = 4, units = "in", res = 300)
hist(wrs.p, main = my_lab,
     xlab = "p value",
     breaks = 50)
dev.off()

