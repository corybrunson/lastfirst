# session library
library(tidyverse)
#remotes::install_github("peekxc/simplextree",
#                        ref = "6e34926b3e6c7c990226e23ba075e29e9a374365")
library(simplextree)
# source simplicial homology code
source(here::here("code/simplicial-homology.r"))
# source directory
rt_data <- "~/Desktop/rt-data" # laptop
#rt_data <- "/blue/rlaubenbacher/jason.brunson/rt-data" # HPG
# store directory
#save_dir <- "/blue/rlaubenbacher/jason.brunson/lastfirst"
save_dir <- "data/cover"

inv_which <- function(i, len) `[<-`(logical(len), i, TRUE)

# function to calculate area under the receiving operator characteristic curve
auc_fun <- function(response, predictor) {
  as.double(suppressMessages(pROC::auc(pROC::roc(
    response = response,
    predictor = predictor
  ))))
}

# extract parameters from file names
list.files(here::here("data/cover/"), full.names = TRUE) %>%
  enframe(name = NULL, value = "file") %>%
  mutate(
    careunit = str_replace(file, "^.*-([A-Z]+)-.*$", "\\1"),
    n_lmks = as.integer(str_replace(file, "^.*-lmk([0-9]+)-.*$", "\\1")),
    ext_mult = as.double(str_replace(file, "^.*-ext([0-9]+)-.*$", "\\1")) / 100,
    proc = str_replace(file, "^.*-([a-z]+)\\.rds$", "\\1")
  ) %>%
  mutate(
    part_coef = NA_real_,
    auc = NA_real_,
    euler = NA_integer_,
    rich = NA_integer_
  ) %>%
  print() -> eval_data

pb <- progress::progress_bar$new(total = nrow(eval_data))
for (i in seq(nrow(eval_data))) {
  pb$tick()
  
  # care unit data
  careunit <- eval_data$careunit[[i]]
  if (i == 1L || careunit != eval_data$careunit[[i - 1L]]) {
    # binary data with cross-validation indices
    file.path(rt_data, str_c("mimic-", tolower(careunit), "-cases.rds")) %>%
      read_rds() ->
      unit_cases
    # binary response matrices
    unit_cases %>%
      select(mortality_hosp) %>%
      mutate_all(as.integer) %>%
      as.matrix() ->
      unit_resp
  }
  # landmark-cover data
  unit_lmks <- read_rds(eval_data$file[[i]])
  
  # partition coefficients
  unit_count <- nrow(unit_cases)
  unit_lmks$cover_set %>%
    # partition matrix
    sapply(inv_which, len = unit_count) %>%
    # weights
    sweep(1L, apply(., 1L, sum), "/") ->
    unit_part
  # partition coefficient
  eval_data$part_coef[[i]] <- sum(unit_part^2) / unit_count
  
  # outcome risk estimates (weighted sum of cover set incidence rates)
  unit_est <- unit_part %*% matrix(unit_lmks$mort)
  # discriminability of risk estimates
  eval_data$auc[[i]] <- auc_fun(as.vector(unit_resp), as.vector(unit_est))
  
  # calculate nerve
  unit_nerve <- nerve(
    simplex_tree(),
    unit_lmks$cover_set,
    k = nrow(unit_lmks) + 1L# k - 1
  )
  # calculate Euler characteristic
  eval_data$euler[[i]] <- euler_characteristic(unit_nerve)
  # calculate homology groups
  unit_hom <- hom_group_ranks_mod2(unit_nerve)
  unit_betti <- unit_hom$homology
  # homological richness
  eval_data$rich[[i]] <- -unit_betti[[1L]] + sum(unit_betti[-1L])
  
  write_rds(eval_data, here::here("data/eval-data.rds"))
}
