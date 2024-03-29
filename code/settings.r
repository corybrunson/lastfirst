# document dimensions (starting with text width obtained using LaTeX)
textwidth <- grid::unit(15.11293 * 4/5, "cm")

# figure dimensions (cm)
textwidth <- 15.11293 * 3/5
phi <- (1 + sqrt(5)) / 2

# global options
options(tibble.print_min = 12L, tibble.print_max = 24L)
#options(ggplot2.discrete.colour = "Set1", ggplot2.discrete.fill = "Set1")
ggplot2::theme_set(ggplot2::theme_bw())
ggplot2::theme_update(legend.position = "bottom", legend.box = "vertical")

# recurring color palettes
proc_pal <- RColorBrewer::brewer.pal(n = 3L, "Set1")
impl_pal <- RColorBrewer::brewer.pal(n = 2L, "Dark2")
frac_pal <- RColorBrewer::brewer.pal(n = 5L, "Accent")[c(seq(3L), 5L)]
# unit_pal <-
#   colorspace::darken(RColorBrewer::brewer.pal(n = 5L, "Pastel2"), .15)
unit_pal <- RColorBrewer::brewer.pal(n = 5L, "Set3")
