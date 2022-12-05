library(tidyverse)

# Higher than average variability position
POS <- 81197167

load("../dataDerived/modelInputData.RData")

samples.df <- read_csv("../dataRaw/masterSamplesheet.csv", show_col_types = F)
covariates.df <- samples.df %>%
  dplyr::mutate(isMale = ifelse(sex == "Male", 1, 0),
                hasLOAD = ifelse(diagnostic_group == "LOAD", 1, 0),
                age = age_at_visit,
                sample_id = as.character(sample_id)) %>%
  dplyr::select(c(sample_id, age, isMale, hasLOAD))



M.df <- M.sub %>%
  dplyr::filter(chromStart == POS) %>%
  pivot_longer(cols = -chromStart, values_to="M", names_to="sample_id") %>%
  dplyr::select(-chromStart)


C.df <- Cov.sub %>%
  dplyr::filter(chromStart == POS) %>%
  pivot_longer(cols = -chromStart, values_to="C", names_to="sample_id") %>%
  dplyr::select(-chromStart)

locus.df <- left_join(M.df, C.df, by = "sample_id") %>%
  left_join(covariates.df, by = "sample_id") %>%
  dplyr::transmute(sample_id, M, C, isMale, hasLOAD,
                   age_std = (age - mean(age)) / sd(age))


save(locus.df, file="../dataDerived/singleLocus.RData")
