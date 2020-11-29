# Python GUDHI setup
library(reticulate)
use_condaenv("r-reticulate")
gd <- import("gudhi")

# convert from simplextree to GUDHI
simplextree_as_gudhi <- function(st) {
  stopifnot(inherits(st, "Rcpp_SimplexTree"))
  lst <- st$as_list()
  gst <- gd$SimplexTree()
  for (i in seq(length(lst))) {
    if (i == 1L) {
      apply(lst[[i]], 2L, function(r) gst$insert(list(r), filtration = i-1L))
    } else {
      apply(lst[[i]], 2L, function(r) gst$insert(r, filtration = i-1L))
    }
  }
  return(gst)
}

# format GUDHI persistence data as a data frame
persistence_to_data_frame <- function(x, pair.col = FALSE) {
  stopifnot(is.list(x))
  res <- data.frame(dimension = vapply(x, `[[`, 1L, i = 1L))
  if (pair.col) {
    res$pair <- lapply(x, `[[`, i = 2L)
  } else {
    res$birth <- vapply(x, `[[`, 1.0, i = c(2L, 1L))
    res$death <- vapply(x, `[[`, 1.0, i = c(2L, 2L))
  }
  res
}
