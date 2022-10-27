# session library
library(tidyverse)

# source and store directories
if (dir.exists("/blue")) {
  # HiPerGator
  mx_data <- "/blue/rlaubenbacher/jason.brunson/covid19-mx/data"
  lastfirst_dir <- "~/lastfirst"
  save_dir <- "/blue/rlaubenbacher/jason.brunson/lastfirst/data/cover"
  # sleep intervals
  sleep_sec <- 0
} else if (str_detect(here::here(), "jason.brunson")) {
  # laptop
  mx_data <- "~/Desktop/covid19-mx/data"
  lastfirst_dir <- here::here()
  save_dir <- "data/cover"
  # sleep intervals
  sleep_sec <- 15
}

# source settings
source(file.path(lastfirst_dir, "code/settings.r"))

# bind all results
file.path(lastfirst_dir, "data") %>%
  list.files("^auc-stats-mx-week-[0-9]+\\.rds", full.names = TRUE) %>%
  enframe(name = NULL, value = "file") %>%
  mutate(data = map(file, read_rds)) %>%
  select(-file) %>%
  unnest(c(data)) %>%
  # mutate(sampler = fct_inorder(sampler)) %>%
  #mutate(sampler = fct_explicit_na(sampler, na_level = "none")) %>%
  mutate(sampler = factor(
    sampler,
    levels = c("maxmin", "lastfirst", "random")
  )) %>%
  mutate(wt_opt = factor(
    wt_opt,
    levels = c("rank1", "rank2", "triangle", "inverse", "gaussian")
  )) %>%
  print() -> auc_stats
# join week cohort sizes
read_rds(file.path(mx_data, "mx-cases.rds")) %>%
  transmute(semana = as.integer(lubridate::week(fecha_ingreso))) %>%
  filter(semana >= min(auc_stats$semana) & semana <= max(auc_stats$semana)) %>%
  group_by(semana) %>% count(name = "size") %>% ungroup() %>%
  right_join(auc_stats, by = "semana") %>%
  print() -> auc_stats
# means
auc_stats %>%
  group_by(semana, size, outcome, sampler, landmarks) %>%
  summarize_at(vars(opt_auc, test_auc), mean) %>%
  ungroup() %>%
  print() -> auc_means

# summarize optimal weighting schemes
auc_stats %>%
  filter(! is.na(sampler)) %>%
  group_by(outcome, sampler, wt_opt) %>%
  count() %>%
  ggplot(aes(x = wt_opt, y = n, fill = sampler)) +
  facet_wrap( ~ outcome) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = proc_pal)

# Gaussian weights only

# bind all results
file.path(lastfirst_dir, "data") %>%
  list.files("^auc-stats-mx-week-gaussian-[0-9]+\\.rds", full.names = TRUE) %>%
  enframe(name = NULL, value = "file") %>%
  mutate(data = map(file, read_rds)) %>%
  select(-file) %>%
  unnest(c(data)) %>%
  # mutate(sampler = fct_inorder(sampler)) %>%
  #mutate(sampler = fct_explicit_na(sampler, na_level = "none")) %>%
  mutate(sampler = factor(
    sampler,
    levels = c("maxmin", "lastfirst", "random")
  )) %>%
  print() -> auc_stats
# join week cohort sizes
read_rds(file.path(mx_data, "mx-cases.rds")) %>%
  transmute(semana = as.integer(lubridate::week(fecha_ingreso))) %>%
  filter(semana >= min(auc_stats$semana) & semana <= max(auc_stats$semana)) %>%
  group_by(semana) %>% count(name = "size") %>% ungroup() %>%
  right_join(auc_stats, by = "semana") %>%
  print() -> auc_stats
# translate and factorize
auc_stats %>%
  mutate(semana = fct_reorder(str_c("Week ", semana), semana)) %>%
  mutate(outcome = fct_recode(
    fct_inorder(outcome),
    Intubation = "intubado", Death = "def"
  )) %>%
  print() -> auc_stats
# means
auc_stats %>%
  group_by(semana, size, outcome, sampler, landmarks) %>%
  summarize_at(vars(opt_auc, test_auc), mean) %>%
  ungroup() %>%
  print() -> auc_means

# summarize optimal settings by sampler, number of landmarks, and week
auc_stats %>%
  #filter(sampler != "none") %>%
  filter(! is.na(sampler)) %>%
  ggplot(aes(x = k_opt)) +
  facet_grid(sampler ~ semana) +
  geom_histogram(binwidth = 12L, position = "dodge") +
  labs(x = "Optimal neighborhood size")

# compare performance across samplers, numbers of landmarks, and week
auc_stats %>%
  ggplot(aes(x = landmarks, y = test_auc, color = sampler)) +
  facet_grid(outcome ~ semana, scales = "fixed") +
  geom_line(data = filter(auc_means, ! is.na(sampler)), size = .5) +
  scale_color_manual(values = proc_pal) +
  labs(x = "Number of landmarks", y = "AUROC", color = "Procedure") ->
  auc_stats_plot
ggsave(here::here("docs/figures/knn-auc-mx-gaussian.pdf"),
       auc_stats_plot,
       width = textwidth, height = textwidth / phi)

# compare performance across samplers, no. landmarks, and week (& w/o sampling)
auc_stats %>%
  ggplot(aes(x = landmarks, y = test_auc, color = sampler)) +
  facet_grid(outcome ~ semana) +
  geom_line(data = filter(auc_means, ! is.na(sampler)), size = .5) +
  geom_hline(data = filter(auc_means, is.na(sampler)),
             aes(yintercept = test_auc), size = .5, color = "darkgrey") +
  scale_color_manual(values = proc_pal) +
  labs(x = "Number of landmarks", y = "AUROC", color = "Procedure") ->
  auc_stats_plot
ggsave(here::here("docs/figures/knn-auc-mx-gaussian.pdf"),
       auc_stats_plot,
       width = textwidth, height = textwidth / phi)
