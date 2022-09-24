x <- tdaunif::sample_circle(n = 24L)

mm <- landmarks_maxmin(
  x, tower = TRUE,
  extend_radius = extension(mult = 1, add = .1)
)
lf <- landmarks_lastfirst(
  x, tower = TRUE,
  extend_cardinality = extension(mult = 1, add = 6L)
)

maxs_to_betti <- function(maxs) {
  st <- simplextree::simplex_tree(maxs)
  stpy <- interplex::as_py_gudhi(st)
  stpy$compute_persistence()
  stpy$betti_numbers()
}
maxs_to_betti(mm$tower_max[[1L]])

mm$tower_max %>%
  sapply(maxs_to_betti) %>%
  t() %>%
  `colnames<-`(c("beta_0", "beta_1")) ->
  mm_betti
lf$tower_max %>%
  sapply(maxs_to_betti) %>%
  t() %>%
  `colnames<-`(c("beta_0", "beta_1")) ->
  lf_betti

table(mm_betti[, "beta_1"])
table(lf_betti[, "beta_1"])
