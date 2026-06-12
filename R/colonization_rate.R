#' Estimates colonization rates events per lineage per unit time
#'
#' For each tip, counts the number of dispersal events on the path from the
#' root to the tip and divides by the total branch length of that path. Each
#' dispersal event is recorded as both an import and an export, so the total
#' transitions divided by two to avoid double-counting. Accounts for incomplete
#' sampling fraction when available.
#'
#' @param phy The marched phylogenetic tree.
#' @param events The data frame returned by \code{map_insitu_events}.
#'
#' @return A data frame with one row per tip and columns:
#'   \describe{
#'     \item{\code{species}}{Tip label.}
#'     \item{\code{n_import}}{Number of import events on the path from root to tip.}
#'     \item{\code{n_export}}{Number of export events on the path from root to tip.}
#'     \item{\code{n_colonization}}{Number of dispersal events}
#'     \item{\code{path_length}}{Total branch length from root to tip.}
#'     \item{\code{colonization_rate}}{Colonization events per unit time.}
#'   }
#' @export
colonization_rate <- function(phy, events) {
  import_nodes <- events$node[!is.na(events$import)]
  export_nodes <- events$node[!is.na(events$export)]
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
    n_import <- sum(path_nodes %in% import_nodes)
    n_export <- sum(path_nodes %in% export_nodes)
    n_colonization <- (n_import + n_export) / 2
    pl <- sum(path_bl) / sampling_fraction
    data.frame(
      species = phy$tip.label[tip_idx],
      n_import = n_import,
      n_export = n_export,
      n_colonization = n_colonization,
      colonization_rate = n_colonization / pl,
      stringsAsFactors  = FALSE
    )
  })

  do.call(rbind, out)
}
