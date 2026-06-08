#' Dated phylogeny with branches and nodes colored by in-situ vs. colonization classification
#'
#' @param phy The phylogenetic tree associated with your data
#' @param events The dataframe returned by `map_insitu_events`. If the user
#' would like to input a custom dataframe, ensure that it has the following
#' columns: node (internal node numbers associated with the given phylogenetic
#' tree), island (the name of the island associated with the ancestral state
#' for that node), and in_situ (logical - whether in situ speciation occurred
#' at the given node).
#' @export

plot_classified_phylo <- function(phy, events){
  islands <- unique(events$island)
  colors  <- setNames(palette.colors(length(islands), "Okabe-Ito"), islands)

  plot(phy,
       main = expression(paste("Island-Level ", italic("in situ"), " Speciation")),
       cex = 0.8)

  # Find the nodes at which in situ speciation occurred
  in_situ_nodes <- events$node[events$in_situ == TRUE]

  # Get the island names for the nodes at which in situ speciation occurred
  node_labels <- events$island[events$node %in% in_situ_nodes]

  # Add points to the nodes where in situ speciation occurred
  ape::nodelabels(node = in_situ_nodes,
                  pch = 19,
                  col  = colors[node_labels])

  # Add island names to the nodes with in situ speciation
  ape::nodelabels(text = node_labels,
                  node = in_situ_nodes,
                  cex = 0.8,
                  frame = "none",
                  adj = c(1.1, -0.4))
}
