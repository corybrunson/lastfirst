# session library
library(tidyverse)
library(binda)

# source and store directories
if (str_detect(here::here(), "corybrunson")) {
  # laptop
  rt_data <- "~/Desktop/rt-data"
  save_dir <- "data/cover"
} else if (str_detect(here::here(), "jason.brunson")) {
  # HiPerGator
  rt_data <- "/blue/rlaubenbacher/jason.brunson/rt-data"
  save_dir <- "/blue/rlaubenbacher/jason.brunson/lastfirst/data/cover"
} else {
  stop("Cannot recognize working directory.")
}

# all care units
careunits <- read_rds(file.path(rt_data, "mimic-units.rds"))
careunits <- c("CCU", "MICU")
# numbers of landmarks
ns_lmks <- c(12L, 60L)
# multiplicative extensions
exts_mult <- 0
# landmark procedures
lmk_procs <- c("mm", "lf")

# minimum BDA score to keep predictor
score_min <- 60
# convenience functions for piping workflows
bda <- function(data, predictors, grouping, ...) {
  stopifnot(
    is.data.frame(data),
    is.character(grouping) || is.numeric(grouping),
    length(grouping) == 1L
  )
  if (is.numeric(grouping)) grouping <- names(data)[[grouping]]
  if (missing(predictors)) predictors <- setdiff(names(data), grouping)
  x <- as.matrix(data[, predictors, drop = FALSE])
  l <- data[[grouping]]
  binda(x, l, ...)
}
bda.rk <- function(data, predictors, grouping, ...) {
  stopifnot(
    is.data.frame(data),
    is.character(grouping) || is.numeric(grouping),
    length(grouping) == 1L
  )
  if (is.numeric(grouping)) grouping <- names(data)[[grouping]]
  if (missing(predictors)) predictors <- setdiff(names(data), grouping)
  x <- as.matrix(data[, predictors, drop = FALSE])
  l <- data[[grouping]]
  binda.ranking(x, l, ...)
}

# compare spreads of mortaliry risks across cover sets

# load all landmark covers
crossing(
  careunit = careunits,
  n_lmks = ns_lmks, lmk_proc = fct_inorder(lmk_procs)
) %>%
  # load landmark covers
  mutate(file = pmap_chr(
    list(careunit, n_lmks, lmk_proc),
    ~ str_c("unit-cover-", ..1, "-lmk", ..2, "-ext0-", ..3, ".rds")
  )) %>%
  mutate(filepath = map_chr(file, ~ file.path(save_dir, .))) %>%
  mutate(data = map(filepath, read_rds)) %>%
  select(-file, -filepath) %>%
  unnest(data) ->
  lmks_data

# compare histograms
lmks_data %>%
  ggplot(aes(x = mort)) +
  facet_grid(lmk_proc ~ careunit + n_lmks) +
  geom_histogram()
# compare standard deviations
lmks_data %>%
  group_by(careunit, n_lmks, lmk_proc) %>%
  summarize(mort_sd = sd(mort)) %>%
  pivot_wider(
    id_cols = c(careunit, n_lmks), names_from = lmk_proc, values_from = mort_sd
  ) %>%
  ungroup()

# binary discriminant analysis (BDA)

# BDA example: CCU with 12 landmarks

# binary data
file.path(rt_data, str_c("mimic-ccu-cases.rds")) %>%
  read_rds() ->
  ccu_cases
# binary predictor and response matrices
ccu_pred <- ccu_cases %>%
  select(-subject_id, -hadm_id, -contains("mortality")) %>%
  mutate_all(as.integer) %>%
  as.matrix()
# column variables
names(ccu_cases) %>%
  enframe(name = "idx", value = "predictor") ->
  ccu_cols

# load unit landmark covers
tibble(lmk_proc = fct_inorder(lmk_procs)) %>%
  # load landmark covers
  mutate(file = map_chr(
    lmk_proc,
    ~ str_c("unit-cover-ccu-lmk12-ext0-", ..1, ".rds")
  )) %>%
  mutate(filepath = map_chr(file, ~ file.path(save_dir, .))) %>%
  mutate(data = map(filepath, read_rds)) %>%
  select(-file, -filepath) %>%
  unnest(data) ->
  lmks_data

# binary discriminant analysis over three highest-mortality cover sets
# (`map()` requires consistent dimensions)
lmks_data %>%
  mutate(label = str_c(lmk_proc, ordinal, sep = "-")) %>%
  group_by(lmk_proc) %>%
  slice_max(order_by = mort, n = 2L) %>%
  ungroup() %>%
  # binary discriminant analysis
  mutate(cover_pred = map(cover_set, ~ ccu_pred[., , drop = FALSE])) %>%
  mutate(cover_pred = map(cover_pred, as_tibble, .name_repair = "unique")) %>%
  select(lmk_proc, label, cover_pred) %>%
  unnest(cover_pred) %>%
  nest(data = -c(lmk_proc)) ->
  ccu_bda_data

# fit
ccu_bda_data %>%
  mutate(cover_bda = map(
    data,
    ~ bda(.x, grouping = "label", verbose = FALSE)
  )) %>%
  select(-data) ->
  ccu_bda_fit
# case probabilities
# -+- TBD -+-

# ranks
ccu_bda_data %>%
  mutate(cover_bda = map(
    data,
    ~ bda.rk(.x, grouping = "label", verbose = FALSE)
  )) %>%
  select(-data) ->
  ccu_bda_rank
# predictor scores
ccu_bda_rank %>%
  mutate(cover_bda = map(cover_bda, ~ as_tibble(unclass(.)))) %>%
  mutate(cover_bda = map(cover_bda, ~ pivot_longer(
    ., cols = starts_with("t."), names_to = "cover_set", values_to = "t.score"
  ))) %>%
  unnest(cover_bda) %>%
  mutate(cover_set = str_remove(cover_set, "^t\\.")) %>%
  mutate(cover_set_idx = as.integer(str_extract(cover_set, "[0-9]+$"))) %>%
  filter(score >= score_min) %>%
  left_join(ccu_cols, by = "idx") ->
  ccu_bda_scores

# comparison of distributions of scores: histogram
ccu_bda_scores %>%
  ggplot(aes(x = score, fill = lmk_proc)) +
  geom_histogram()
# comparison of distributions of scores: size-rank plot
ccu_bda_scores %>%
  group_by(lmk_proc) %>%
  arrange(score) %>%
  mutate(rank = seq(1L, length(score))) %>%
  ggplot(aes(x = rank, y = score, color = lmk_proc)) +
  geom_line()

# comparison of specific predictors
ccu_bda_scores %>%
  filter(lmk_proc == "mm") %>%
  mutate(predictor = fct_reorder(predictor, score)) %>%
  ggplot(aes(x = t.score, y = predictor, fill = cover_set)) +
  geom_col() +
  guides(fill = "none") ->
  ccu_12_mm_bda_scores
ccu_bda_scores %>%
  filter(lmk_proc == "lf") %>%
  mutate(predictor = fct_reorder(predictor, score)) %>%
  ggplot(aes(x = t.score, y = predictor, fill = cover_set)) +
  geom_col() +
  guides(fill = "none") ->
  ccu_12_lf_bda_scores
ccu_bda_scores %>%
  mutate(predictor = str_c(ifelse(lmk_proc == "mm", "", " "), predictor)) %>%
  mutate(cover_set = fct_reorder(cover_set, as.integer(lmk_proc))) %>%
  ggplot(aes(x = t.score, y = fct_reorder(predictor, score), fill = cover_set)) +
  facet_grid(lmk_proc ~ ., scales = "free_y") +
  geom_col() +
  labs(x = "t-score", y = "Predictor", fill = "Cover set")

# BDA across all settings

# initialize data set
bda_scores <- tibble()
# loop over care units
for (careunit in careunits) {
  
  # binary data
  file.path(rt_data, str_c("mimic-", tolower(careunit), "-cases.rds")) %>%
    read_rds() ->
    unit_cases
  # binary predictor and response matrices
  unit_pred <- unit_cases %>%
    select(-subject_id, -hadm_id, -contains("mortality")) %>%
    mutate_all(as.integer) %>%
    as.matrix()
  # column variables
  names(unit_cases) %>%
    enframe(name = "idx", value = "predictor") ->
    unit_cols
  
  # load unit landmark covers
  crossing(n_lmks = ns_lmks, lmk_proc = fct_inorder(lmk_procs)) %>%
    # load landmark covers
    mutate(file = pmap_chr(
      list(n_lmks, lmk_proc),
      ~ str_c("unit-cover-", careunit, "-lmk", ..1, "-ext0-", ..2, ".rds")
    )) %>%
    mutate(filepath = map_chr(file, ~ file.path(save_dir, .))) %>%
    mutate(data = map(filepath, read_rds)) %>%
    select(-file, -filepath) %>%
    unnest(data) ->
    lmks_data
  
  # binary discriminant analysis over three highest-mortality cover sets
  # (`map()` requires consistent dimensions)
  lmks_data %>%
    mutate(label = str_c(careunit, n_lmks, lmk_proc, ordinal, sep = "-")) %>%
    group_by(n_lmks, lmk_proc) %>%
    slice_max(order_by = mort, n = 2L) %>%
    ungroup() %>%
    # binary discriminant analysis
    mutate(cover_pred = map(cover_set, ~ unit_pred[., , drop = FALSE])) %>%
    mutate(cover_pred = map(cover_pred, as_tibble, .name_repair = "unique")) %>%
    select(n_lmks, lmk_proc, label, cover_pred) %>%
    unnest(cover_pred) %>%
    nest(data = -c(n_lmks, lmk_proc)) ->
    unit_bda_data
  
  # fit
  unit_bda_data %>%
    mutate(cover_bda = map(
      data,
      ~ bda(.x, grouping = "label", verbose = FALSE)
    )) %>%
    select(-data) ->
    unit_bda_fit
  # case probabilities
  
  
  # ranks
  unit_bda_data %>%
    mutate(cover_bda = map(
      data,
      ~ bda.rk(.x, grouping = "label", verbose = FALSE)
    )) %>%
    select(-data) ->
    unit_bda_rank
  # predictor scores
  unit_bda_rank %>%
    mutate(cover_bda = map(cover_bda, ~ as_tibble(unclass(.)))) %>%
    mutate(cover_bda = map(cover_bda, ~ pivot_longer(
      ., cols = starts_with("t."), names_to = "cover_set", values_to = "t.score"
    ))) %>%
    unnest(cover_bda) %>%
    mutate(cover_set = str_remove(cover_set, "^t\\.")) %>%
    mutate(cover_set_idx = as.integer(str_extract(cover_set, "[0-9]+$"))) %>%
    filter(score >= score_min) %>%
    left_join(unit_cols, by = "idx") ->
    unit_bda_scores
  bda_scores <- bind_rows(
    bda_scores,
    mutate(unit_bda_scores, careunit = careunit)
  )
}
