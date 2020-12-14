# session library
library(tidyverse)
library(landmark)

# source and store directories
if (str_detect(here::here(), "corybrunson")) {
  # laptop
  rt_data <- "~/Desktop/rt-data"
  lastfirst_dir <- here::here()
  save_dir <- "data/cover"
  # sleep intervals
  sleep_sec <- 15
} else if (str_detect(here::here(), "jason.brunson")) {
  # HiPerGator
  rt_data <- "/blue/rlaubenbacher/jason.brunson/rt-data"
  lastfirst_dir <- "~/lastfirst"
  save_dir <- "/blue/rlaubenbacher/jason.brunson/lastfirst/data/cover"
  # sleep intervals
  sleep_sec <- 0
} else {
  stop("Cannot recognize working directory.")
}

# folds in cross-validation
o_folds <- 6L
i_folds <- 6L
# range of neighborhood size to consider
max_k <- 180L
min_k <- 12L
# weight functions take a k-by-test matrix and apply a function to each column
# (the matrix contains cosine similarity values ranging from -1 to 1)
# https://github.com/KlausVigo/kknn/blob/master/R/kknn.R
wt_funs <- list(
  rank1 = function(d) (nrow(d) + 1L) - apply(d, 2L, rank, ties.method = "min"),
  rank2 = function(d) (nrow(d) + 1L) - apply(d, 2L, rank, ties.method = "max"),
  triangle = function(d) 1 - d,
  inverse = function(d) 1 / d,
  gaussian = function(d) {
    alpha <- 1 / (2 * (nrow(d) + 1))
    qua <- abs(qnorm(alpha))
    t(apply(d * qua, 1, dnorm, sd = 1))
  }
)
# numbers of landmarks
ns_lmks <- c(36L, 60L, 180L, 360L)
# landmark generators
lmk_funs <- list(
  random = function(x, num) {
    set.seed(7392L) # randomly-generated seed
    sample(nrow(x), num)
  },
  maxmin = function(x, num) {
    landmarks_maxmin(x, num = num,
                     dist_method = "cosine", seed_index = "minmax")
  },
  lastfirst = function(x, num) {
    landmarks_lastfirst(x, num = num,
                        dist_method = "cosine", seed_index = "firstlast")
  }
)
# function to calculate area under the receiving operator characteristic curve
auc_fun <- function(response, predictor) {
  as.double(suppressMessages(pROC::auc(pROC::roc(
    response = response,
    predictor = predictor
  ))))
}

if (file.exists(file.path(lastfirst_dir, "data/auc-stats-micu.rds"))) {
  readr::read_rds(file.path(lastfirst_dir, "data/auc-stats-micu.rds")) ->
    auc_stats
} else {
  # initialize data frame
  auc_stats <- tibble()
}

# add stratified cross-validation indices to binary data
file.path(rt_data, str_c("mimic-micu-cases.rds")) %>%
  read_rds() %>%
  group_by(mortality_hosp) %>%
  mutate(row = row_number()) %>%
  mutate(outer = (sample(row) %% o_folds) + 1L) %>%
  group_by(mortality_hosp, outer) %>%
  mutate(row = row_number()) %>%
  mutate(inner = (sample(row) %% i_folds) + 1L) %>%
  select(-row) %>%
  ungroup() %>%
  print() -> unit_cases
unit_cases %>% select(mortality_hosp, outer, inner) %>% table() %>% print()

# loop over folds
for (i in seq(o_folds)) for (j in seq(i_folds)) {
  n_ij <- length(which(auc_stats$outer == i & auc_stats$inner == j))
  if (n_ij < length(ns_lmks) * length(lmk_funs)) {
    # begin this incomplete round from scratch
    auc_stats <- filter(auc_stats, outer != i | inner != j)
  } else {
    # skip this complete round
    next
  }
  
  # training, optimizing, and testing indices
  train <- which(unit_cases$outer != i & unit_cases$inner != j)
  opt <- which(unit_cases$outer != i & unit_cases$inner == j)
  test <- which(unit_cases$outer == i)
  
  # binary predictor and response matrices
  unit_pred <- unit_cases %>%
    select(-subject_id, -hadm_id, -outer, -inner,
           -contains("mortality")) %>%
    mutate_all(as.integer) %>%
    as.matrix()
  unit_resp <- unit_cases %>%
    select(mortality_hosp) %>%
    mutate_all(as.integer) %>%
    as.matrix()
  
  print(c(i, j))
  
  # nearest training neighbors of each optimizing datum
  nbrs <- proxy::dist(unit_pred[train, ], unit_pred[opt, ],
                      method = "cosine") %>%
    unclass() %>% as.data.frame() %>% as.list() %>% unname() %>%
    lapply(enframe, name = "id", value = "dist") %>%
    lapply(mutate, rank = rank(dist, ties.method = "min")) %>%
    lapply(filter, rank <= max_k)
  # predictions for each landmark using each neighborhood size
  k_opt_preds <- sapply(nbrs, function(df) {
    vapply(seq(max_k),
           function(k) mean(unit_resp[train[df$id[df$rank <= k]], ]),
           FUN.VALUE = 1)
  })
  
  # predictive accuracy on optimizing data
  k_aucs <- apply(k_opt_preds, 1L, auc_fun, response = unit_resp[opt, ])
  
  # (first) best parameters subject to neighborhoods of size at least 12
  max_auc <- max(k_aucs[seq(min_k, length(k_aucs))])
  k_opt <- which(k_aucs[seq(min_k, length(k_aucs))] == max_auc,
                 arr.ind = TRUE)[1L] + (min_k - 1L)
  
  # predictions on testing data
  test_preds <- proxy::dist(unit_pred[train, ], unit_pred[test, ],
                            method = "cosine") %>%
    unclass() %>% as.data.frame() %>% as.list() %>% unname() %>%
    lapply(enframe, name = "id", value = "dist") %>%
    lapply(mutate, rank = rank(dist, ties.method = "min")) %>%
    lapply(filter, rank <= k_opt) %>%
    sapply(function(df) mean(unit_resp[train[df$id[df$rank <= k_opt]], ]))
  test_auc <- auc_fun(response = unit_resp[test, ],
                      predictor = as.vector(test_preds))
  
  # augment data frame
  auc_stats <- bind_rows(auc_stats, tibble(
    careunit = "MICU",
    outer = i, inner = j,
    sampler = NA_character_, landmarks = NA_integer_,
    k_wt_auc = list(k_aucs),
    k_opt = k_opt, wt_opt = NA_character_,
    opt_auc = max_auc, test_auc = test_auc
  ))
  Sys.sleep(sleep_sec)
  
  # loop over landmark-generating functions and numbers of landmarks
  for (l in seq_along(lmk_funs)) for (n_lmks in ns_lmks) {
    print(c(i, j, names(lmk_funs)[[l]], n_lmks))
    
    # training set landmarks
    lmks <- lmk_funs[[l]](unit_pred[train, ], num = n_lmks)
    # nearest training neighbors of these landmarks
    nbrs <- proxy::dist(unit_pred[train, ], unit_pred[train[lmks], ],
                        method = "cosine") %>%
      unclass() %>% as.data.frame() %>% as.list() %>% unname() %>%
      lapply(enframe, name = "id", value = "dist") %>%
      lapply(mutate, rank = rank(dist, ties.method = "min")) %>%
      lapply(filter, rank <= max_k)
    # predictions for each landmark using each neighborhood size
    k_lmk_preds <- sapply(nbrs, function(df) {
      vapply(seq(max_k),
             function(k) mean(unit_resp[train[df$id[df$rank <= k]], ]),
             FUN.VALUE = 1)
    })
    
    # predictive accuracy on optimizing data
    lmk_opt_dists <- proxy::dist(unit_pred[train[lmks], ], unit_pred[opt, ],
                                 method = "cosine")
    k_wt_aucs <- array(NA_real_, dim = c(max_k, length(wt_funs)))
    k_wt_auc_data <- k_wt_aucs %>%
      as_tibble() %>%
      set_names(names(wt_funs)) %>%
      mutate(k = row_number())
    
    # loop over weighting functions
    for (w in seq_along(wt_funs)) {
      opt_k_preds <- t(k_lmk_preds %*% wt_funs[[w]](lmk_opt_dists)) /
        colSums(wt_funs[[w]](lmk_opt_dists))
      opt_auc <- apply(opt_k_preds, 2L,
                       auc_fun, response = unit_resp[opt, ])
      k_wt_aucs[, w] <- opt_auc
    }
    
    # (first) best parameters subject to neighborhoods of size at least 12
    max_auc <- max(k_wt_aucs[seq(min_k, nrow(k_wt_aucs)), ])
    k_wt_opt <- which(k_wt_aucs[seq(min_k, nrow(k_wt_aucs)), ] == max_auc,
                      arr.ind = TRUE)[1L, ] + c(min_k - 1L, 0L)
    
    # predictions on testing data
    lmk_test_dists <- proxy::dist(unit_pred[train[lmks], ],
                                  unit_pred[test, ],
                                  method = "cosine")
    test_preds <- k_lmk_preds[k_wt_opt[[1L]], , drop = FALSE] %*%
      wt_funs[[k_wt_opt[[2L]]]](lmk_test_dists) /
      colSums(wt_funs[[k_wt_opt[[2L]]]](lmk_test_dists))
    test_auc <- auc_fun(response = unit_resp[test, ],
                        predictor = as.vector(test_preds))
    
    # augment data frame
    auc_stats <- bind_rows(auc_stats, tibble(
      careunit = "MICU",
      outer = i, inner = j,
      sampler = names(lmk_funs)[[l]], landmarks = n_lmks,
      k_wt_auc = list(k_wt_auc_data),
      k_opt = k_wt_opt[[1L]], wt_opt = names(wt_funs)[[k_wt_opt[[2L]]]],
      opt_auc = max_auc, test_auc = test_auc
    ))
    Sys.sleep(sleep_sec)
    
  }
  
  write_rds(auc_stats, file.path(lastfirst_dir, "data/auc-stats-micu.rds"))
}
