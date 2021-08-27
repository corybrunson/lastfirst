# session library
library(tidyverse)
library(simplextree)

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

inv_which <- function(i, len) `[<-`(logical(len), i, TRUE)

# function to calculate area under the receiving operator characteristic curve
auc_fun <- function(response, predictor) {
  as.double(suppressMessages(pROC::auc(pROC::roc(
    response = response,
    predictor = predictor
  ))))
}

# Euler characteristic
euler_characteristic <- function(st) {
  sum(st$n_simplices * rep_len(c(1, -1), st$dimension + 1L))
}

# all careunits
careunits <- read_rds(file.path(rt_data, "mimic-units.rds"))

if (file.exists(file.path(lastfirst_dir, "data/eval-data.rds"))) {
  eval_data <- read_rds(file.path(lastfirst_dir, "data/eval-data.rds"))
} else {
  # extract parameters from file names
  list.files(save_dir, full.names = FALSE) %>%
    enframe(name = NULL, value = "file") %>%
    filter(str_detect(file, "^unit-cover")) %>%
    mutate(
      careunit = toupper(str_replace(file, "^.*-([a-zA-Z]+)-.*$", "\\1")),
      n_lmks = as.integer(str_replace(file, "^.*-lmk([0-9]+)-.*$", "\\1")),
      ext_mult = as.double(str_replace(file, "^.*-ext([0-9]+)-.*$", "\\1")) /
        100,
      proc = str_replace(file, "^.*-([a-z]+)\\.rds$", "\\1")
    ) %>%
    mutate(
      careunit = factor(careunit, levels = intersect(careunits, careunit)),
      proc = factor(proc, levels = c("mm", "lf"))
    ) %>%
    mutate(
      part_coef = NA_real_,
      auc = NA_real_,
      simplices = vector("list", nrow(.)),
      euler = NA_integer_,
      #rich = NA_integer_,
      done = FALSE
    ) %>%
    arrange(careunit, n_lmks, ext_mult, proc) ->
    eval_data
}

undone_rows <- which(! eval_data$done)
pb <- progress::progress_bar$new(total = length(undone_rows))
for (i in undone_rows) {
  message("Evaluating results in row ", i, ":")
  pb$tick()
  
  message("Reading care unit data...")
  # care unit data
  careunit <- eval_data$careunit[[i]]
  if (! exists("unit_cases") ||
      i == 1L || careunit != eval_data$careunit[[i - 1L]]) {
    # binary response matrix
    file.path(rt_data, str_c("mimic-", tolower(careunit), "-cases.rds")) %>%
      read_rds() %>%
      select(mortality_hosp) %>%
      mutate_all(as.integer) %>%
      as.matrix() ->
      unit_resp
    unit_count <- nrow(unit_resp)
  }
  
  message("Reading landmarks data...")
  # landmark-cover data
  unit_lmks <- read_rds(file.path(save_dir, eval_data$file[[i]]))
  
  message("Computing partition coefficient...")
  # partition coefficient
  unit_lmks$cover_set %>%
    # partition matrix
    sapply(inv_which, len = unit_count) %>%
    # weights
    sweep(1L, apply(., 1L, sum), "/") ->
    unit_part
  # partition coefficient
  eval_data$part_coef[[i]] <- sum(unit_part^2) / unit_count
  
  message("Computing discriminability statistic...")
  # outcome risk estimates (weighted sum of cover set incidence rates)
  unit_est <- unit_part %*% matrix(unit_lmks$mort)
  # discriminability of risk estimates
  eval_data$auc[[i]] <- auc_fun(as.vector(unit_resp), as.vector(unit_est))
  rm(unit_part)
  
  message("Computing cover nerve...")
  # calculate nerve
  unit_nerve <- nerve(
    simplex_tree(),
    unit_lmks$cover_set,
    k = 1L + 1L# k - 1
  )
  # count simplices of each dimension
  eval_data$simplices[[i]] <- unit_nerve$n_simplices
  # calculate Euler characteristic
  eval_data$euler[[i]] <- euler_characteristic(unit_nerve)
  # calculate homology groups
  #unit_hom <- hom_group_ranks_mod2(unit_nerve)
  #unit_betti <- unit_hom$homology
  # homological richness
  #eval_data$rich[[i]] <- -unit_betti[[1L]] + sum(unit_betti[-1L])
  
  eval_data$done[[i]] <- TRUE
  message("Writing evaluation data...")
  write_rds(eval_data, file.path(lastfirst_dir, "data/eval-data.rds"))
}
