#' Dated phylogeny with branches and nodes colored by in-situ vs. colonization classification
#'
#' @param phy The phylogenetic tree associated with your data
#' @param events The dataframe returned by `map_insitu_events`. If the user
#' would like to input a custom dataframe, ensure that it has the following
#' columns: node (internal node numbers associated with the given phylogenetic
#' tree), island (the name of the island associated with the ancestral state
#' for that node), and in_situ (logical - whether in situ speciation occurred
#' at the given node).
#' @param PAM The presence-absence matrix associated with your data as returned
#' by the `match_island_phylo` function.
#' @export

plot_classified_phylo <- function(phy, events, PAM){
  islands <- PAM$locale
  colors <- stats::setNames(grDevices::palette.colors(length(islands), "Okabe-Ito"), islands)

  # Plot the basic phylogeny, but add a label offset for the pie charts and
  #  some extra space for the legend (based on tree depth to make it general)
  plot(phy,
       main = expression(paste("Island-Level ", italic("in situ"), " Speciation")),
       label.offset = 0.5,
       x.lim = c(0, max(ape::node.depth.edgelength(phy)) * 1.5),
       cex = 0.8)

  # Find the nodes at which in situ speciation occurred
  in_situ_nodes <- events$node[events$in_situ == TRUE]

  # Get the island names for the nodes at which in situ speciation occurred
  node_labels <- events$island[events$node %in% in_situ_nodes]

  # Add points to the nodes where in situ speciation occurred
  ape::nodelabels(node = in_situ_nodes,
                  pch = 19,
                  col = colors[node_labels])

  # Add island names to the nodes with in situ speciation
  ape::nodelabels(text = node_labels,
                  node = in_situ_nodes,
                  cex = 0.8,
                  frame = "none",
                  adj = c(1.1, -0.4))

  # Add pie charts to tip labels to show extant islands
  # Make sure that the locales are the row names so that they don't get removed
  #  when the columns are reordered
  rownames(PAM) <- PAM$locale

  # Reorder the PAM columns to match tip labels
  # Inspiration: https://stackoverflow.com/questions/25446714/reorder-matrix-columns-by-matching-column-names-to-vector-of-names
  PAM <- PAM[, phy$tip.label, drop = FALSE]

  # The tiplabels() function needs species names to be row names
  #  to add pie charts
  pie_mat <- t(PAM)

  # Get proportions of island occurrence for each species
  pie_mat <- pie_mat / rowSums(pie_mat)

  # Add pie charts to tips
  ape::tiplabels(pie = pie_mat,
                 piecol = colors[colnames(pie_mat)],
                 cex = 0.3)

  # Add legend, with a location based on the tree size in both dimensions
  graphics::legend(x = max(ape::node.depth.edgelength(phy)) * 1.3,
                   y = ape::Ntip(phy) / 2,
                   legend = islands,
                   pch = 19,
                   col = colors,
                   bty = "n",
                   cex = 0.7)
}
