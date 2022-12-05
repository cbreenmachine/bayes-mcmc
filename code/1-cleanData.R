library(tidyverse)
library(data.table)

#CEP131 gene
#chr17:81,189,593-81,222,999
left <- 81189593
right <- 81222999

# Read in data
M <- fread("../dataRaw/chr17.M.bed")
Cov <- fread("../dataRaw/chr17.Cov.bed")

# Filter to only include CpGs with minimum 5x coverage
vv <- apply(X=M[ , -"chromStart"] / Cov[ , -"chromStart"], FUN=var, MARGIN=1)
cc <- apply(X=Cov[ , -"chromStart"], FUN=min, MARGIN=1)
keepix <- (cc >= 5 & vv > 0.01)
sum(keepix)


df <- read_csv("../dataRaw/masterSamplesheet.csv", show_col_types = F)

ten_loci <- M$chromStart[keepix][1:10]

# Subset to one region
# Pick the first ten CpGs
M.sub <- M %>% dplyr::filter(chromStart %in% ten_loci)
Cov.sub <- Cov %>% dplyr::filter(chromStart %in% ten_loci)

# Check order
all(names(M.sub) == names(Cov.sub))

X <- df %>%
  dplyr::filter(sample_id %in% colnames(M.sub)) %>%
  distinct() %>%
  column_to_rownames("sample_id") %>%
  dplyr::transmute(
    load_status = ifelse(diagnostic_group == "LOAD", 1, 0),
    age_std = (age_at_visit - mean(age_at_visit)) / sd(age_at_visit),
    is_male = ifelse(sex == "Male", 1, 0))

# WIll be used in RSTan
y <- dplyr::select(M.sub, rownames(X))
C <- dplyr::select(Cov.sub, rownames(X))

# Check one last time that everything matches
all(rownames(X) == colnames(y))

K <- nrow(y)
P <- ncol(X)
N <- nrow(X)

save(list = c("y", "C", "X", "K", "P", "N"), file = "../dataDerived/20221130-tenLoci-v1.RData")

