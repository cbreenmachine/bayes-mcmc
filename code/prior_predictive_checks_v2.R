library(tidyverse)
library(ggsci)

load("../dataDerived/20221130-tenLoci-v1.RData")


ebayes.mu <- mean(unlist(y / C))
ebayes.nu <- ebayes.mu * (1 - ebayes.mu) / 2

all_coverage <- unlist(C)
