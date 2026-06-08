#' Plot island node classification on a phylogeny
#'
#' @param phy The matched phylogenetic tree .
#' @param node_class The data frame returned by \code{identify_island_nodes}.
#' @param PAM A presence-absence matrix.
#' @export
plot_island_nodes <- function(phy, node_class, PAM) {
  islands <- unique(PAM$locale)
  colors  <- setNames(palette.colors(length(islands), "Okabe-Ito"), islands)

  plot(phy,
       main = "Island node classification",
       cex  = 0.8)

  # within-island nodes colored by island
  within <- node_class[node_class$status == "within_island", ]
  ape::nodelabels(node = within$node, pch = 19, col = colors[within$island])

  # between-island nodes in grey
  between <- node_class[node_class$status == "between_island", ]
  ape::nodelabels(node = between$node, pch = 19, col = "grey60")

  legend("bottomleft",
         legend = c(islands, "between-island"),
         pch    = 19,
         col    = c(colors, "grey60"),
         bty    = "n",
         cex    = 0.7)
}
