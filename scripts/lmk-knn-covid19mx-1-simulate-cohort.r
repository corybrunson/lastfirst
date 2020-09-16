
# session library
library(tidyverse)
library(landmark)
# source directories
covid19mx_data <- "~/Desktop/covid19-mx/data"

# load pre-processed data
file.path("~/Desktop/covid19-mx/data/covid19mx.rds") %>%
  read_rds() %>%
  print() -> covid19mx_cases

# sleep intervals
sleep_sec <- 15
# number of patients per cohort
cohort_size <- 1000L
# folds in cross-validation
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

# subset of data for predictive modeling
covid19mx_cases %>%
  filter(resultado != "resultado pendiente") %>%
  filter(fecha_ingreso > as.Date("2020-03-10")) %>%
  arrange(fecha_ingreso) %>%
  print() -> covid19mx_subset

# split for consecutive-window CV (at least 1000 patients per cohort)
covid19mx_subset %>%
  group_by(entidad_um, fecha_ingreso) %>%
  count() %>%
  group_by(entidad_um) %>%
  mutate(n_cohorts = as.integer(floor(sum(n) / cohort_size))) %>%
  mutate(cohort = as.integer(ceiling(cumsum(n) / cohort_size))) %>%
  filter(cohort <= n_cohorts) %>%
  select(entidad_um, fecha_ingreso, cohort) %>%
  ungroup() %>%
  print() -> covid19mx_cohorts

# process data for modeling
covid19mx_subset %>%
  # flag and restrict to cohorts
  inner_join(covid19mx_cohorts, by = c("entidad_um", "fecha_ingreso")) %>%
  # split for within-window optimization
  group_by(entidad_um, fecha_ingreso, resultado) %>%
  mutate(row = row_number()) %>%
  mutate(inner = (sample(row) %% i_folds) + 1L) %>%
  select(-row) %>%
  ungroup() %>%
  # un-factor `entidad_um`
  mutate(entidad_um = as.character(entidad_um)) %>%
  # collapse infrequent sectors to 'other'
  mutate(sector = fct_lump_n(sector, n = 60L, other_level = "otra")) %>%
  # two age bracket variables, to increase influence and reduce arbitrariness
  mutate(edad_a_0 = cut(edad, breaks = seq(-6L, max(edad) + 10L, 10L)),
         edad_a_5 = cut(edad, breaks = seq(-1L, max(edad) + 10L, 10L))) %>%
  select(-edad) %>%
  mutate_if(is.factor, fct_drop) %>%
  fastDummies::dummy_cols(
    c("sector", "sexo", "entidad_res",
      "edad_a_5", "edad_a_0",
      "tipo_paciente", "nacionalidad"),
    remove_selected_columns = TRUE
  ) %>%
  mutate_at(vars(intubado, neumonia, diabetes:otro_caso, uci),
            ~ ifelse(. == "si", 1L, 0L)) %>%
  mutate(resultado = ifelse(resultado == "positivo sars-cov-2", 1L, 0L)) %>%
  print() -> covid19mx_model
# only `entidad_um`, dates, and durations are non-integers
print(which(! sapply(covid19mx_model, is.integer)))
# `cohort` and `inner` are only non-predictor integers
print(which(sapply(select(covid19mx_model, -entidad_um), max) != 1L))

# initialize list and data frame
auc_stats <- tibble()

# loop over states
for (entidad in unique(covid19mx_model$entidad_um)) {
  
  # loop over folds
  for (i in unique(covid19mx_model$cohort)[-1L]) for (j in seq(i_folds)) {
    
    # training, optimizing, and testing indices
    train <- which(covid19mx_model$cohort == i - 1L)
    opt <- which(covid19mx_model$cohort == i &
                   covid19mx_model$inner == j)
    test <- which(covid19mx_model$cohort == i &
                    covid19mx_model$inner != j)
    
    # binary predictor and response matrices
    covid19mx_pred <- covid19mx_model %>%
      select_if(is.integer) %>%
      select(-cohort, -inner, -resultado) %>%
      as.matrix()
    covid19mx_resp <- covid19mx_model %>%
      select(resultado) %>%
      as.matrix()
    
    print(c(entidad, i, j))
    
    # nearest training neighbors of each optimizing datum
    nbrs <- proxy::dist(covid19mx_pred[train, ], covid19mx_pred[opt, ],
                        method = "cosine") %>%
      unclass() %>% as.data.frame() %>% as.list() %>% unname() %>%
      lapply(enframe, name = "id", value = "dist") %>%
      lapply(mutate, rank = rank(dist, ties.method = "min")) %>%
      lapply(filter, rank <= max_k)
    # predictions for each landmark using each neighborhood size
    k_opt_preds <- sapply(nbrs, function(df) {
      vapply(seq(max_k),
             function(k) mean(covid19mx_resp[train[df$id[df$rank <= k]], ]),
             FUN.VALUE = 1)
    })
    
    # predictive accuracy on optimizing data
    k_aucs <- apply(k_opt_preds, 1L, auc_fun, response = covid19mx_resp[opt, ])
    
    # (first) best parameters subject to neighborhoods of size at least 12
    max_auc <- max(k_aucs[seq(min_k, length(k_aucs))])
    k_opt <- which(k_aucs[seq(min_k, length(k_aucs))] == max_auc,
                   arr.ind = TRUE)[1L] + (min_k - 1L)
    
    # predictions on testing data
    test_preds <- proxy::dist(covid19mx_pred[train, ], covid19mx_pred[test, ],
                              method = "cosine") %>%
      unclass() %>% as.data.frame() %>% as.list() %>% unname() %>%
      lapply(enframe, name = "id", value = "dist") %>%
      lapply(mutate, rank = rank(dist, ties.method = "min")) %>%
      lapply(filter, rank <= k_opt) %>%
      sapply(function(df)
        mean(covid19mx_resp[train[df$id[df$rank <= k_opt]], ]))
    test_auc <- auc_fun(response = covid19mx_resp[test, ],
                        predictor = as.vector(test_preds))
    
    # augment data frame
    auc_stats <- bind_rows(auc_stats, tibble(
      entidad = entidad, cohort = i, inner = j,
      sampler = NA_character_, landmarks = NA_integer_,
      k_wt_auc = list(k_aucs),
      k_opt = k_opt, wt_opt = NA_character_,
      opt_auc = max_auc, test_auc = test_auc
    ))
    Sys.sleep(sleep_sec)
    
    # loop over landmark-generating functions and numbers of landmarks
    for (l in seq_along(lmk_funs)) for (n_lmks in ns_lmks) {
      
      print(c(entidad, i, j, names(lmk_funs)[[l]], n_lmks))
      
      # training set landmarks
      lmks <- lmk_funs[[l]](covid19mx_pred[train, ], num = n_lmks)
      # nearest training neighbors of these landmarks
      nbrs <- proxy::dist(covid19mx_pred[train, ],
                          covid19mx_pred[train[lmks], ],
                          method = "cosine") %>%
        unclass() %>% as.data.frame() %>% as.list() %>% unname() %>%
        lapply(enframe, name = "id", value = "dist") %>%
        lapply(mutate, rank = rank(dist, ties.method = "min")) %>%
        lapply(filter, rank <= max_k)
      # predictions for each landmark using each neighborhood size
      k_lmk_preds <- sapply(nbrs, function(df) {
        vapply(seq(max_k),
               function(k) mean(covid19mx_resp[train[df$id[df$rank <= k]], ]),
               FUN.VALUE = 1)
      })
      
      # predictive accuracy on optimizing data
      lmk_opt_dists <- proxy::dist(covid19mx_pred[train[lmks], ],
                                   covid19mx_pred[opt, ],
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
                         auc_fun, response = covid19mx_resp[opt, ])
        k_wt_aucs[, w] <- opt_auc
      }
      
      # (first) best parameters subject to neighborhoods of size at least 12
      max_auc <- max(k_wt_aucs[seq(min_k, nrow(k_wt_aucs)), ])
      k_wt_opt <- which(k_wt_aucs[seq(min_k, nrow(k_wt_aucs)), ] == max_auc,
                        arr.ind = TRUE)[1L, ] + c(min_k - 1L, 0L)
      
      # predictions on testing data
      lmk_test_dists <- proxy::dist(covid19mx_pred[train[lmks], ],
                                    covid19mx_pred[test, ],
                                    method = "cosine")
      test_preds <- k_lmk_preds[k_wt_opt[[1L]], , drop = FALSE] %*%
        wt_funs[[k_wt_opt[[2L]]]](lmk_test_dists) /
        colSums(wt_funs[[k_wt_opt[[2L]]]](lmk_test_dists))
      test_auc <- auc_fun(response = covid19mx_resp[test, ],
                          predictor = as.vector(test_preds))
      
      # augment data frame
      auc_stats <- bind_rows(auc_stats, tibble(
        entidad = entidad, cohort = i, inner = j,
        sampler = names(lmk_funs)[[l]], landmarks = n_lmks,
        k_wt_auc = list(k_wt_auc_data),
        k_opt = k_wt_opt[[1L]], wt_opt = names(wt_funs)[[k_wt_opt[[2L]]]],
        opt_auc = max_auc, test_auc = test_auc
      ))
      Sys.sleep(sleep_sec)
      
    }
    
  }
  
  write_rds(auc_stats, "data/auc-stats-covid19mx-cohort.rds")
  
}
