#' In-situ speciation index per island
#'
#' Computes a standardized index of in-situ contribution to diversity for
#' each island. The index (isci) is the number of in-situ speciation events
#' on an island divided by the total branch length of the island's subtree.
#' This is, the sum of all branch lengths of all descendant tips that are
#' exclusively assigned to that island.
#'
#' @param phy Matched phylogenetic tree
#' @param events The data frame returned by \code{map_insitu_events}
#'
#' @return A data frame with one row per island and columns:
#'   \describe{
#'     \item{\code{island}}{Island name.}
#'     \item{\code{n_insitu}}{Number of in-situ speciation events.}
#'     \item{\code{island_bl}}{Total branch length of the island subtree.}
#'     \item{\code{isci}}{In-situ contribution index: \code{n_insitu / island_bl}.}
#'   }
#' @export

insitu_speciation_index <- function(phy, events){
  events <- events[events$island != "TRANSITION",]
  islands <- unique(events$island[!is.na(events$island)])
  n_tips  <- ape::Ntip(phy)

  out <- lapply(islands, function(isl) {
    # nodes assigned to this island
    isl_nodes <- events$node[events$island == isl]

    # branch lengths of edges whose child node is in isl_nodes
    edge_mask <- phy$edge[, 2] %in% isl_nodes
    island_bl <- sum(phy$edge.length[edge_mask])

    # in-situ events on this island
    n_insitu <- sum(events$in_situ[events$island == isl], na.rm = TRUE)

    data.frame(
      island    = isl,
      n_insitu  = n_insitu,
      island_bl = island_bl,
      isci      = if (island_bl > 0) n_insitu / island_bl else NA,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}
