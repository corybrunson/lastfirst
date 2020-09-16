
# session library
library(tidyverse)

# bind all results
here::here("data/auc-stats.rds") %>%
  read_rds() %>%
  mutate(sampler = fct_explicit_na(sampler, na_level = "none")) %>%
  print() -> auc_stats

# compare performance across samplers, numbers of landmarks, and care units
auc_stats %>%
  filter(sampler != "none") %>%
  ggplot(aes(x = factor(landmarks), y = test_auc, color = sampler)) +
  coord_flip() +
  facet_wrap(~ careunit) +
  geom_boxplot()

# compare performance to basic nearest-neighbors prediction
auc_stats %>%
  ggplot(aes(x = factor(landmarks), y = test_auc, color = sampler)) +
  coord_flip() +
  facet_wrap(~ careunit) +
  geom_boxplot()

# compare optimal neighborhood size across samplers, landmarks, and care units
auc_stats %>%
  filter(sampler != "none") %>%
  ggplot(aes(x = factor(landmarks), y = k_opt, color = sampler)) +
  coord_flip() +
  facet_wrap(~ careunit) +
  geom_boxplot()
