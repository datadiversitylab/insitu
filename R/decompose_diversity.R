#' Partition island species richness into in-situ, import, and export components
#'
#' For each island, counts in-situ speciation, import, and export events.
#' In-situ speciation events are identified as cases in which an ancestral
#' node is assigned to the same region as its descendant. Import events are
#' cases where a descendant node is assigned to a different region than its
#' ancestor. Export events are the reverse: the ancestral node is assigned
#' to the focal island while its descendant is assigned to a different region.
#' Both ancestor and descendant must be assigned to a region for any event
#' to be detected.
#'
#' @param phy The matched phylogenetic tree.
#' @param events The data frame returned by \code{map_insitu_events}, with
#'   columns \code{node}, \code{island}, \code{in_situ}, \code{export}, and
#'   \code{import}.
#' @param PAM The presence-absence matrix with islands as rows and species as
#'   columns, with a \code{locale} column for island names.
#'
#' @return A data frame with one row per island and columns:
#'   \describe{
#'     \item{\code{island}}{Island name.}
#'     \item{\code{n_species}}{Total extant species on the island.}
#'     \item{\code{n_insitu}}{Species derived from in-situ speciation.}
#'     \item{\code{n_import}}{Import events into this island.}
#'     \item{\code{n_export}}{Export events out of this island.}
#'     \item{\code{prop_insitu}}{Proportion of species from in-situ speciation.}
#'   }
#' @export
decompose_diversity <- function(phy, events, PAM) {
  species_cols <- setdiff(names(PAM), "locale")

  tips_from_node <- function(nd) {
    desc  <- integer(0)
    queue <- nd
    while (length(queue) > 0) {
      children <- phy$edge[phy$edge[, 1] == queue[1], 2]
      tips <- children[children <= ape::Ntip(phy)]
      internal <- children[children >  ape::Ntip(phy)]
      desc <- c(desc, tips)
      queue <- c(queue[-1], internal)
    }
    phy$tip.label[desc]
  }

  islands <- unique(PAM$locale)

  out <- lapply(islands, function(isl) {
    isl_tips <- species_cols[sapply(species_cols, function(sp) PAM[[sp]][PAM$locale == isl] == 1)]
    isl_tips <- isl_tips[isl_tips %in% phy$tip.label]
    n_species <- length(isl_tips)

    if (n_species == 0L) {
      return(data.frame(island = isl, n_species = 0L, n_insitu = 0L,
                        n_import = 0L, n_export = 0L, prop_insitu = NA_real_,
                        stringsAsFactors = FALSE))
    }

    # n_insitu: unique tips descending from in-situ nodes on this island
    insitu_nodes <- events$node[events$island == isl & events$in_situ == TRUE]
    insitu_tips <- unique(unlist(lapply(insitu_nodes, tips_from_node)))
    n_insitu <- length(intersect(insitu_tips, isl_tips))

    # n_import: events where this island is the destination
    n_import <- sum(events$import == isl, na.rm = TRUE)

    # n_export: events where this island is the source
    n_export <- sum(events$export == isl, na.rm = TRUE)

    data.frame(
      island = isl,
      n_species = n_species,
      n_insitu = n_insitu,
      n_import = n_import,
      n_export = n_export,
      prop_insitu = n_insitu / n_species,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}
