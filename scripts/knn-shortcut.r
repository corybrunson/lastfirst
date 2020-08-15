
# session library
library(tidyverse)
library(landmark)
# source directories
rt_data <- "~/Desktop/rt-data"
# all care units
careunits <- read_rds(file.path(rt_data, "mimic-units.rds"))
careunits <- setdiff(careunits, c("NICU", "NWARD"))[1]

# folds in cross-validation
o_folds <- 6L
i_folds <- 6L
# range of neighborhood size to consider
max_k <- 180L
min_k <- 12L
# weight functions take a matrix and apply a function to each column
# https://github.com/KlausVigo/kknn/blob/master/R/kknn.R
wt_funs <- list(
  rank1 = function(d) (nrow(d) + 1L) - apply(d, 2L, rank, ties.method = "min"),
  rank2 = function(d) (nrow(d) + 1L) - apply(d, 2L, rank, ties.method = "max"),
  inverse = function(d) 1 / d,
  triangle = function(d) 1 - d
)
# numbers of landmarks
ns_lmks <- c(30L, 60L, 120L)
auc_fun <- function(response, predictor) {
  as.double(suppressMessages(pROC::auc(pROC::roc(
    response = response,
    predictor = predictor
  ))))
}

# uniform 6-by-6-fold cross-validation split
folds <- data.frame(outer = rep(seq(o_folds), times = i_folds),
                    inner = rep(seq(i_folds), each = o_folds))

# initialize list and data frame
res_dims <- c(length(careunits), length(ns_lmks), o_folds, i_folds)
auc_plots <- vector(mode = "list",
                    length = prod(res_dims))
dim(auc_plots) <- res_dims
auc_stats <- tibble()

# loop over all care units
for (careunit in careunits) {
  
  # binary data with cross-validation indices
  file.path(rt_data, str_c("mimic-", tolower(careunit), "-cases.rds")) %>%
    read_rds() %>%
    bind_cols(folds[sample(nrow(folds), nrow(.), replace = TRUE), ]) %>%
    print() -> unit_cases
  
  # loop over numbers of landmarks
  for (n_lmks in ns_lmks) {
    
    # loop over folds
    for (i in seq(o_folds)) for (j in seq(i_folds)) {
      print(c(careunit, n_lmks, i, j))
      
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
      
      # training set landmarks
      lmks <- landmarks_maxmin(unit_pred[train, ], dist_method = "cosine",
                               num = n_lmks, seed_index = "minmax")
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
      for (w in seq_along(wt_funs)) {
        opt_k_preds <- t(k_lmk_preds %*% wt_funs[[w]](lmk_opt_dists)) /
          colSums(wt_funs[[w]](lmk_opt_dists))
        opt_auc <- apply(opt_k_preds, 2L, auc_fun, response = unit_resp[opt, ])
        k_wt_aucs[, w] <- opt_auc
      }
      
      # plot AUC across neighborhood size by weight function!
      auc_plot <- k_wt_aucs %>%
        as_tibble() %>%
        set_names(names(wt_funs)) %>%
        mutate(k = row_number()) %>%
        gather(-k, key = "Weight", value = "AUC") %>%
        ggplot(aes(x = k, y = AUC, group = Weight)) +
        geom_line(aes(color = Weight)) +
        ggtitle(str_c("Optimization step (", i, ",", j, ") for ",
                      careunit, " with ", n_lmks, " landmarks"))
      auc_plots[[match(careunit, careunits), match(n_lmks, ns_lmks), i, j]] <-
        auc_plot
      # (first) best parameters subject to neighborhoods of size at least 12
      max_auc <- max(k_wt_aucs[seq(min_k, nrow(k_wt_aucs)), ])
      k_wt_opt <- which(k_wt_aucs[seq(min_k, nrow(k_wt_aucs)), ] == max_auc,
                        arr.ind = TRUE)[1, ] + c(min_k - 1L, 0L)
      
      # predictions on testing data
      lmk_test_dists <- proxy::dist(unit_pred[train[lmks], ], unit_pred[test, ],
                                    method = "cosine")
      test_preds <- k_lmk_preds[k_wt_opt[[1L]], , drop = FALSE] %*%
        wt_funs[[k_wt_opt[[2L]]]](lmk_test_dists) /
        colSums(wt_funs[[k_wt_opt[[2L]]]](lmk_test_dists))
      test_auc <- auc_fun(response = unit_resp[test, ],
                          predictor = as.vector(test_preds))
      # augment data frame
      auc_stats <- bind_rows(auc_stats, tibble(
        careunit = careunit, landmarks = n_lmks,
        outer = i, inner = j,
        k_opt = k_wt_opt[[1L]], wt_opt = names(wt_funs)[[k_wt_opt[[2L]]]],
        opt_auc = max_auc, test_auc = test_auc
      ))
    }
    
  }
}

write_rds(auc_stats, "data/auc-stats.rds")
write_rds(auc_plots, "data/auc-plots.rds")
