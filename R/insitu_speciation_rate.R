#' Estimates the rate of in-situ speciation events per lineage per unit time
#'
#' For each tip, counts the number of in-situ speciation events on the path
#' from the root to the tip and divides that by the total branch length of the
#' path. When the tree was prepared with \code{\link{prep_phylo}} and a
#' \code{sampling_fraction} less than 1 was set, branch lengths are scaled
#' upward by dividing by the sampling fraction before computing rates.
#' This approach corrects for the fact that incomplete sampling shortens
#' the total branch length of the observed tree relative to the true underlying tree.
#'
#' @param phy The matched phylogenetic tree prepared with \code{\link{prep_phylo}}
#' @param events The data frame returned by \code{map_insitu_events}.
#'
#' @return A data frame with one row per tip and columns:
#'   \describe{
#'     \item{\code{species}}{Tip label.}
#'     \item{\code{n_insitu}}{Number of in-situ speciation events on the path from root to tip.}
#'     \item{\code{path_length}}{Total branch length from root to tip.}
#'     \item{\code{insitu_rate}}{In-situ speciation events per unit time.}
#'   }
#' @export
insitu_speciation_rate <- function(phy, events) {
  insitu_nodes      <- events$node[events$in_situ == TRUE]
  sampling_fraction <- attr(phy, "sampling_fraction")
  if (is.null(sampling_fraction)) sampling_fraction <- 1

  out <- lapply(seq_len(ape::Ntip(phy)), function(tip_idx) {
    path_nodes <- integer(0)
    path_bl <- numeric(0)
    current <- tip_idx
    while (TRUE) {
      parent_row <- which(phy$edge[, 2] == current)
      if (length(parent_row) == 0L) break
      path_nodes <- c(path_nodes, current)
      path_bl <- c(path_bl, phy$edge.length[parent_row])
      current <- phy$edge[parent_row, 1]
    }
    n_insitu <- sum(path_nodes %in% insitu_nodes)
    pl <- sum(path_bl) / sampling_fraction
    data.frame(
      species          = phy$tip.label[tip_idx],
      n_insitu         = n_insitu,
      insitu_rate      = n_insitu / pl,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}
