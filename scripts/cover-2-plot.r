# session library
library(tidyverse)
#remotes::install_github("peekxc/simplextree",
#                        ref = "6e34926b3e6c7c990226e23ba075e29e9a374365")
library(simplextree)
# source directory
rt_data <- "~/Desktop/rt-data" # laptop
#rt_data <- "/blue/rlaubenbacher/jason.brunson/rt-data" # HPG
# store directory
#save_dir <- "/blue/rlaubenbacher/jason.brunson/lastfirst"
save_dir <- "data/cover"

# calculate simplex summaries from case values
summarize_nerve <- function(cover, values, funs, k = 1L) {
  # construct nerve from cover
  nrv <- nerve(simplex_tree(), cover, k = k)
  # calculate simplex dimensions
  n_k_simplices <- nrv$n_simplices[seq(k + 1L)]
  nrv_dimensions <- rep(seq(0L, k), n_k_simplices)
  # index of each simplex (within dimension)
  nrv_ids <- unlist(lapply(nrv$n_simplices, seq)) - 1L
  # traversals for the simplex dimensions
  nrv_traversals <- lapply(seq(0L, k), k_simplices, st = nrv)
  # list observations in each simplex
  nrv_observations <- unlist(lapply(nrv_traversals, ltraverse, function(sigma) {
    Reduce(intersect, cover[sigma])
  }), recursive = FALSE)
  # count observations in each simplex
  nrv_sizes <- vapply(nrv_observations, length, 0L)
  # data frame of basic summary stats
  summ1 <- data.frame(
    dimension = nrv_dimensions,
    ordinal = nrv_ids,
    size = nrv_sizes
  )
  # ensure that `funs` is a list
  if (! is.list(funs)) {
    stopifnot(is.function(funs))
    funs <- list(funs)
  } else {
    stopifnot(all(sapply(funs, is.function)))
  }
  # data frame of additional summary functions
  summ2 <- lapply(funs, function(f) {
    sapply(lapply(nrv_observations, function(x) values[x]), f)
  })
  names(summ2) <- names(funs)
  summ2 <- as.data.frame(summ2)
  # bind data frames
  cbind(summ1, summ2)
}

# all care units
careunits <- read_rds(file.path(rt_data, "mimic-units.rds"))
careunits <- setdiff(careunits, c("NICU", "NWARD"))
# numbers of landmarks
ns_lmks <- c(12L, 24L, 36L, 48L, 60L, 120L, 180L)
# multiplicative extensions
exts_mult <- c(0, .1, .25)
# landmark procedures
lmk_procs <- c("mm", "lf")

careunit <- careunits[[1L]]
n_lmks <- ns_lmks[[1L]]
ext_mult <- exts_mult[[1L]]
lmk_proc <- lmk_procs[[1L]]

# load landmark cover
unit_lmks <- read_rds(file.path(save_dir, str_c(
  "unit-cover-", careunit, "-lmk", n_lmks, "-ext", ext_mult*100, "-", lmk_proc,
  ".rds"
)))

# calculate nerve (1-skeleton)
(unit_nerve <- nerve(simplex_tree(), unit_lmks$cover_set, k = 1L))
# calculate layout
unit_nerve$as_edge_list() %>%
  igraph::graph_from_edgelist(directed = FALSE) %>%
  igraph::layout_with_fr() ->
  unit_coords
# summarize nerve
unit_summ <- summarize_nerve(
  unit_lmks$cover_set, as.vector(unit_resp), list(mort = mean), k = 1L
)
unit_mort <- (unit_summ$mort - min(unit_summ$mort)) /
  (max(unit_summ$mort) - min(unit_summ$mort))
unit_mort <- (unit_lmks$mort - min(unit_lmks$mort)) /
  (max(unit_lmks$mort) - min(unit_lmks$mort))
unit_pal <- rainbow(n = 256L)[1 + unit_mort * 255]
unit_size_0 <- unit_lmks$size
unit_size_1 <- apply(unit_nerve$edges, 1L, function(ij) {
  length(intersect(
    unit_lmks$cover_set[[ij[[1L]]]],
    unit_lmks$cover_set[[ij[[2L]]]]
  ))
})
# plot
plot(
  unit_nerve,
  coords = unit_coords,
  vertex_opt = list(
    pch = 16,
    cex = 6 * unit_size_1 / max(unit_size_1),
    col = unit_pal
  ),
  edge_opt = list(
    lwd = 5 * unit_size_1 / max(unit_size_1),
    col = "darkgrey"
  )
)
