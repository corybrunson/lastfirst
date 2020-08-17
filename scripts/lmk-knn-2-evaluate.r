
# session library
library(tidyverse)

auc_stats <- read_rds("data/auc-stats.rds")

auc_stats %>%
  ggplot(aes(x = factor(landmarks), y = test_auc, group = sampler)) +
  facet_wrap(~ careunit) +
  geom_boxplot()


stop()
# plot AUC across neighborhood size by weight function!
auc_plot <- k_wt_auc_data %>%
  gather(-k, key = "Weight", value = "AUC") %>%
  ggplot(aes(x = k, y = AUC, group = Weight)) +
  geom_line(aes(color = Weight)) +
  ggtitle(str_c("CV (", i, ",", j, ") for ", careunit, " with ",
                n_lmks, " ", names(lmk_funs)[[l]], " landmarks"))
auc_plots[[match(careunit, careunits), match(n_lmks, ns_lmks), i, j]] <-
  auc_plot
