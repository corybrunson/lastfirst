# session library
library(tidyverse)
library(simplextree)

# source and store directories
if (str_detect(here::here(), "corybrunson")) {
  # laptop
  mx_data <- "~/Desktop/covid19-mx/data"
  save_dir <- "data/cover"
  lastfirst_dir <- here::here()
} else if (str_detect(here::here(), "jason.brunson")) {
  # HiPerGator
  mx_data <- "/blue/rlaubenbacher/jason.brunson/covid19-mx/data"
  save_dir <- "/blue/rlaubenbacher/jason.brunson/lastfirst/data/cover"
  lastfirst_dir <- "~/lastfirst"
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

if (file.exists(file.path(lastfirst_dir, "data/eval-data-mx.rds"))) {
  eval_data <- read_rds(file.path(lastfirst_dir, "data/eval-data-mx.rds"))
} else {
  # extract parameters from file names
  list.files(save_dir, full.names = FALSE) %>%
    enframe(name = NULL, value = "file") %>%
    filter(str_detect(file, "^mx-cover")) %>%
    mutate(
      n_lmks = as.integer(str_replace(file, "^.*-lmk([0-9]+)-.*$", "\\1")),
      ext_mult = as.double(str_replace(file, "^.*-ext([0-9]+)-.*$", "\\1")) /
        100,
      proc = str_replace(file, "^.*-([a-z]+)\\.rds$", "\\1")
    ) %>%
    mutate(
      proc = factor(proc, levels = c("mm", "lf"))
    ) %>%
    mutate(
      part_coef = NA_real_,
      auc_intubado = NA_real_, auc_def = NA_real_,
      simplices = vector("list", nrow(.)),
      euler = NA_integer_,
      #rich = NA_integer_,
      done = FALSE
    ) %>%
    arrange(n_lmks, ext_mult, proc) ->
    eval_data
}

# covid-19 data binary response matrix
file.path(mx_data, "mx-cases.rds") %>%
  read_rds() %>%
  select(intubado, def) %>%
  mutate_all(as.integer) %>%
  as.matrix() ->
  mx_resp
mx_count <- nrow(mx_resp)

undone_rows <- which(! eval_data$done)
pb <- progress::progress_bar$new(total = length(undone_rows))
for (i in undone_rows) {
  message("Evaluating results in row ", i, ":")
  pb$tick()
  
  message("Reading landmarks data...")
  # landmark-cover data
  mx_lmks <- read_rds(file.path(save_dir, eval_data$file[[i]]))
  
  message("Computing partition coefficient...")
  # partition coefficient
  mx_lmks$cover_set %>%
    # partition matrix
    sapply(inv_which, len = mx_count) %>%
    # weights
    sweep(1L, apply(., 1L, sum), "/") ->
    mx_part
  # partition coefficient
  eval_data$part_coef[[i]] <- sum(mx_part^2) / mx_count
  
  message("Computing discriminability statistics...")
  # outcome risk estimates (weighted sum of cover set incidence rates)
  mx_intubado_est <- mx_part %*% matrix(mx_lmks$intubado)
  mx_def_est <- mx_part %*% matrix(mx_lmks$def)
  # discriminability of risk estimates
  eval_data$auc_intubado[[i]] <-
    auc_fun(as.vector(mx_resp[, 1L]), as.vector(mx_intubado_est))
  eval_data$auc_def[[i]] <-
    auc_fun(as.vector(mx_resp[, 2L]), as.vector(mx_def_est))
  rm(mx_part)
  
  message("Computing cover nerve...")
  # calculate nerve
  mx_nerve <- nerve(
    simplex_tree(),
    mx_lmks$cover_set,
    k = 1L + 1L# k - 1
  )
  # count simplices of each dimension
  eval_data$simplices[[i]] <- mx_nerve$n_simplices
  # calculate Euler characteristic
  eval_data$euler[[i]] <- euler_characteristic(mx_nerve)
  # calculate homology groups
  #mx_hom <- hom_group_ranks_mod2(mx_nerve)
  #mx_betti <- mx_hom$homology
  # homological richness
  #eval_data$rich[[i]] <- -mx_betti[[1L]] + sum(mx_betti[-1L])
  
  eval_data$done[[i]] <- TRUE
  message("Writing evaluation data...")
  write_rds(eval_data, file.path(lastfirst_dir, "data/eval-data-mx.rds"))
}
