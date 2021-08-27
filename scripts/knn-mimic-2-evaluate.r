# session library
library(tidyverse)

# source and store directories
if (stringr::str_detect(here::here(), "corybrunson")) {
  # laptop
  machine <- "Cory's MacBook Air"
  rt_data <- "~/Desktop/rt-data"
  lastfirst_dir <- here::here()
  library(landmark)
} else if (stringr::str_detect(here::here(), "Users/jason.brunson")) {
  # desktop
  machine <- "Cory's UF iMac"
  rt_data <- "~/Desktop/rt-data"
  lastfirst_dir <- here::here()
  devtools::load_all("~/Documents/proj-active/tda/landmark/")
} else if (stringr::str_detect(here::here(), "home/jason.brunson")) {
  # HiPerGator
  machine <- "HiPerGator cluster"
  rt_data <- "/blue/rlaubenbacher/jason.brunson/rt-data"
  lastfirst_dir <- "~/lastfirst"
  library(landmark)
} else {
  stop("Cannot recognize working directory.")
}

# source settings
source(file.path(lastfirst_dir, "code/settings.r"))

# bind all results
file.path(lastfirst_dir, "data") %>%
  list.files("^auc-stats-[a-z]+\\.rds", full.names = TRUE) %>%
  enframe(name = NULL, value = "file") %>%
  mutate(data = map(file, read_rds)) %>%
  select(-file) %>%
  unnest(c(data)) %>%
  mutate(sampler = fct_inorder(sampler)) %>%
  #mutate(sampler = fct_explicit_na(sampler, na_level = "none")) %>%
  mutate(wt_opt = factor(
    wt_opt,
    levels = c("rank1", "rank2", "triangle", "inverse", "gaussian")
  )) %>%
  print() -> auc_stats
# join care unit sizes
file.path(rt_data) %>%
  list.files("^mimic-[a-zA-Z]+-cases.rds", full.names = TRUE) %>%
  enframe(name = NULL, value = "file") %>%
  mutate(data = map(file, read_rds)) %>%
  transmute(
    careunit = str_replace(file, "^.*mimic-([a-zA-Z]+)-cases.*$", "\\1"),
    careunit = toupper(careunit),
    size = map_int(data, nrow)
  ) %>%
  right_join(auc_stats, by = "careunit") %>%
  # order care units by size
  mutate(careunit = fct_reorder(careunit, size)) %>%
  print() -> auc_stats

# summarize optimal settings by sampler, number of landmarks, and care unit
auc_stats %>%
  #filter(sampler != "none") %>%
  filter(! is.na(sampler)) %>%
  ggplot(aes(x = k_opt)) +
  facet_grid(sampler ~ careunit) +
  geom_histogram(binwidth = 6L, position = "dodge") +
  labs(x = "Optimal neighborhood size")

# compare performance across samplers, numbers of landmarks, and care units
auc_stats %>%
  #filter(sampler != "none") %>%
  filter(! is.na(sampler)) %>%
  ggplot(aes(x = factor(landmarks), y = test_auc, color = sampler)) +
  #coord_flip() +
  facet_wrap(~ careunit) +
  geom_boxplot() +
  labs(x = "Number of landmarks", y = "AUROC", color = "Procedure") ->
  auc_stats_plot
ggsave(here::here("docs/figures/knn-auc-1.pdf"),
       auc_stats_plot,
       width = textwidth, height = textwidth / phi, units = "cm")

# compare performance to basic nearest-neighbors prediction
auc_stats %>%
  ggplot(aes(x = factor(landmarks), y = test_auc, color = sampler)) +
  #coord_flip() +
  facet_wrap(~ careunit) +
  geom_boxplot() +
  labs(x = "Number of landmarks", y = "AUROC", color = "Procedure") ->
  auc_stats_plot
ggsave(here::here("docs/figures/knn-auc-2.pdf"),
       auc_stats_plot,
       width = textwidth, height = textwidth / phi, units = "cm")

auc_stats %>%
  mutate(careunit = str_c(
    careunit, " (n = ", format(size, big.mark = ","), ")", sep = ""
  )) %>%
  mutate(careunit = fct_reorder(careunit, size)) %>%
  mutate(sampler = fct_recode(
    sampler,
    Random = "random",
    `Similarity threshold` = "maxmin",
    `Case count` = "lastfirst"
  )) %>%
  ggplot(aes(x = factor(landmarks), y = test_auc, color = sampler)) +
  facet_grid(~ careunit) +
  geom_boxplot() +
  labs(x = "Number of landmarks",
       y = "AUROC (mortality)",
       color = "Procedure") ->
  auc_stats_plot
ggsave(here::here("docs/figures/knn-auc.jpg"),
       auc_stats_plot,
       width = textwidth, height = textwidth / 4, units = "cm")
