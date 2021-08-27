# session library
library(tidyverse)
library(landmark)

# source and store directories
if (str_detect(here::here(), "corybrunson")) {
  # laptop
  mx_data <- "~/Desktop/covid19-mx/data"
  lastfirst_dir <- here::here()
} else if (str_detect(here::here(), "jason.brunson")) {
  # HiPerGator
  mx_data <- "/blue/rlaubenbacher/jason.brunson/covid19-mx/data"
  lastfirst_dir <- "~/lastfirst"
} else {
  stop("Cannot recognize working directory.")
}

# folds in cross-validation
i_folds <- 6L
# range of neighborhood size to consider
max_k <- 180L
min_k <- 12L
# weight function
gaussian_weight <- function(d) {
  alpha <- 1 / (2 * (nrow(d) + 1))
  qua <- abs(qnorm(alpha))
  t(apply(d * qua, 1, dnorm, sd = 1))
}
# numbers of landmarks
ns_lmks <- c(36L, 48L, 60L, 120L, 180L, 240L)
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

# initialize stats data frame
auc_stats <- tibble()

# Mexican covid-19 data
read_rds(file.path(mx_data, "mx-cases.rds")) %>%
  mutate(semana_ingreso = as.integer(lubridate::week(fecha_ingreso))) %>%
  select(-starts_with("fecha_")) ->
  mx_cases
# semanas
semanas_ingreso <- unique(mx_cases$semana_ingreso)
# outcomes
outcomes <- c("intubado", "def")

for (o in seq_along(outcomes)) {
  outcome <- outcomes[[o]]
  
  # Mexican covid-19 data (overwrite only 'inner' variable each loop)
  mx_cases %>%
    # split for within-window optimization
    group_by_at(vars(all_of(c("semana_ingreso", outcome)))) %>%
    mutate(row = row_number()) %>%
    mutate(inner = (sample(row) %% i_folds) + 1L) %>%
    select(-row) %>%
    ungroup() ->
    mx_cases
  
  # loop over folds
  for (i in semanas_ingreso[-1L]) for (j in seq(i_folds)) {
    message("Fitting for semana ", i, " of ", max(semanas_ingreso),
            " and fold ", j, " of ", i_folds, "...")
    
    # training, optimizing, and testing indices
    train <- which(mx_cases$semana_ingreso == i - 1L)
    opt <- which(mx_cases$semana_ingreso == i & mx_cases$inner == j)
    test <- which(mx_cases$semana_ingreso == i & mx_cases$inner != j)
    
    # binary predictor and response matrices
    mx_pred <- mx_cases %>%
      select(-semana_ingreso, -inner, -intubado, -def) %>%
      mutate_all(as.integer) %>%
      as.matrix()
    mx_resp <- mx_cases %>%
      select(intubado, def) %>%
      mutate_all(as.integer) %>%
      as.matrix()
    
    # nearest training neighbors of each optimizing datum
    nbrs <- proxy::dist(mx_pred[train, ], mx_pred[opt, ],
                        method = "cosine") %>%
      unclass() %>% as.data.frame() %>% as.list() %>% unname() %>%
      lapply(enframe, name = "id", value = "dist") %>%
      lapply(mutate, rank = rank(dist, ties.method = "min")) %>%
      lapply(filter, rank <= max_k)
    # predictions for each landmark using each neighborhood size
    k_opt_preds <- sapply(nbrs, function(df) {
      vapply(seq(max_k),
             function(k) mean(mx_resp[train[df$id[df$rank <= k]], o]),
             FUN.VALUE = 1)
    })
    
    # predictive accuracy on optimizing data
    k_aucs <- apply(k_opt_preds, 1L, auc_fun, response = mx_resp[opt, o])
    
    # (first) best parameters subject to neighborhoods of size at least 12
    max_auc <- max(k_aucs[seq(min_k, length(k_aucs))])
    k_opt <- which(k_aucs[seq(min_k, length(k_aucs))] == max_auc,
                   arr.ind = TRUE)[1L] + (min_k - 1L)
    
    # predictions on testing data
    test_preds <- proxy::dist(mx_pred[train, ], mx_pred[test, ],
                              method = "cosine") %>%
      unclass() %>% as.data.frame() %>% as.list() %>% unname() %>%
      lapply(enframe, name = "id", value = "dist") %>%
      lapply(mutate, rank = rank(dist, ties.method = "min")) %>%
      lapply(filter, rank <= k_opt) %>%
      sapply(function(df) mean(mx_resp[train[df$id[df$rank <= k_opt]], ]))
    test_auc <- auc_fun(response = mx_resp[test, o],
                        predictor = as.vector(test_preds))
    
    # augment data frame
    auc_stats <- bind_rows(auc_stats, tibble(
      outcome = outcome, semana = i, inner = j,
      sampler = NA_character_, landmarks = NA_integer_,
      k_aucs = list(k_aucs), k_opt = k_opt,
      opt_auc = max_auc, test_auc = test_auc
    ))
    
    # loop over landmark-generating functions and numbers of landmarks
    for (l in seq_along(lmk_funs)) for (n_lmks in ns_lmks) {
      message("Fitting with ", names(lmk_funs)[[l]], " selection of ",
              n_lmks, " landmarks...")
      
      # training set landmarks
      lmks <- lmk_funs[[l]](mx_pred[train, ], num = n_lmks)
      # nearest training neighbors of these landmarks
      nbrs <- proxy::dist(mx_pred[train, ], mx_pred[train[lmks], ],
                          method = "cosine") %>%
        unclass() %>% as.data.frame() %>% as.list() %>% unname() %>%
        lapply(enframe, name = "id", value = "dist") %>%
        lapply(mutate, rank = rank(dist, ties.method = "min")) %>%
        lapply(filter, rank <= max_k)
      # predictions for each landmark using each neighborhood size
      k_lmk_preds <- sapply(nbrs, function(df) {
        vapply(seq(max_k),
               function(k) mean(mx_resp[train[df$id[df$rank <= k]], o]),
               FUN.VALUE = 1)
      })
      
      # predictive accuracy on optimizing data
      lmk_opt_dists <- proxy::dist(mx_pred[train[lmks], ],
                                   mx_pred[opt, ],
                                   method = "cosine")
      
      opt_k_preds <- t(k_lmk_preds %*% gaussian_weight(lmk_opt_dists)) /
        colSums(gaussian_weight(lmk_opt_dists))
      k_aucs <- apply(opt_k_preds, 2L, auc_fun, response = mx_resp[opt, o])
      
      # (first) best parameters subject to neighborhoods of size at least 12
      max_auc <- max(k_aucs[seq(min_k, length(k_aucs))])
      k_opt <- which(k_aucs[seq(min_k, length(k_aucs))] == max_auc,
                        arr.ind = TRUE)[1L] + (min_k - 1L)
      
      # predictions on testing data
      lmk_test_dists <- proxy::dist(mx_pred[train[lmks], ],
                                    mx_pred[test, ],
                                    method = "cosine")
      test_preds <- k_lmk_preds[k_opt, , drop = FALSE] %*%
        gaussian_weight(lmk_test_dists) /
        colSums(gaussian_weight(lmk_test_dists))
      test_auc <- auc_fun(response = mx_resp[test, o],
                          predictor = as.vector(test_preds))
      
      # augment data frame
      auc_stats <- bind_rows(auc_stats, tibble(
        outcome = outcome, semana = i, inner = j,
        sampler = names(lmk_funs)[[l]], landmarks = n_lmks,
        k_aucs = list(k_aucs), k_opt = k_opt,
        opt_auc = max_auc, test_auc = test_auc
      ))
      
    }
    
    write_rds(auc_stats, file.path(
      lastfirst_dir,
      str_c("data/auc-stats-mx-week-gaussian-", i, ".rds")
    ))
  }
  
}
