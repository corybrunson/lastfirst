# session library
library(tidyverse)
library(simplextree)


ccu_cover <- read_rds("data/visualization/unit-cover-CCU-lmk12-ext0-lf.rds")
ccu_nerve <- read_rds("data/visualization/unit-nerve-CCU-lmk12-ext0-lf.rds")
#simplex_tree(ccu_nerve)
ccu_st <- simplex_tree()
for (s in ccu_nerve) ccu_st$insert(s)


pdf(file = "test.pdf")
plot(ccu_st)
dev.off()
