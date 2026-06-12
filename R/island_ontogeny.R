#' Flags and filters in-situ events relative to island formation age
#'
#' Compares the age of each classified node to the formation age of the
#' island it is assigned to. Nodes that predate island formation are
#' biologically impossible and represent either errors in the ancestral state
#' reconstruction or lineages whose common ancestor actually predates the island.
#' These nodes are flagged and optionally removed. This function also computes
#' island-age-corrected rates by replacing the full path length denominator
#' with the time elapsed since island formation, which is the relevant
#' timescale for assessing how rapidly an island has diversified.
#'
#' @param events The data frame returned by \code{map_insitu_events}.
#' @param phy The matched phylogenetic tree.
#' @param island_ages A named numeric vector of island formation ages, in
#'   the same units as the branch lengths of \code{phy}. Names must match
#'   the island names in \code{events}.
#' @param remove_predating Logical. If \code{TRUE}, nodes that predate
#'   island formation are removed from the returned events table. If
#'   \code{FALSE} (default), they are retained but flagged in the
#'   \code{predates_island} column.
#'
#' @return The events data frame with three additional columns:
#'   \describe{
#'     \item{\code{node_age}}{Age of the node in the same units as branch
#'       lengths (time before present).}
#'     \item{\code{island_age}}{Formation age of the island assigned to
#'       that node.}
#'     \item{\code{predates_island}}{Logical. \code{TRUE} if the node age
#'       exceeds the island formation age, meaning the event predates the
#'       island's existence.}
#'   }
#'
#' @export
age_correct_events <- function(events, phy, island_ages,
                                remove_predating = FALSE) {
  if (!is.numeric(island_ages) || is.null(names(island_ages)))
    stop("'island_ages' must be a named numeric vector.")

  n_tips <- ape::Ntip(phy)
  depths <- ape::node.depth.edgelength(phy)
  root_age <- max(depths[seq_len(n_tips)])
  node_ages <- root_age - depths

  events$node_age <- node_ages[events$node]
  events$island_age <- island_ages[events$island]
  events$predates_island <- !is.na(events$island_age) &
    events$node_age > events$island_age

  n_flagged <- sum(events$predates_island, na.rm = TRUE)
  if (n_flagged > 0)
    message(n_flagged, " node(s) predate their island's formation age and ",
            if (remove_predating) "have been removed." else "are flagged.")

  if (remove_predating)
    events <- events[!events$predates_island, ]

  events
}


#' Computes island-age-corrected in-situ speciation rates
#'
#' Modifies \code{\link{insitu_speciation_rate}} to use island formation age
#' as the time denominator rather than total root-to-tip path length. This
#' gives the rate of in-situ speciation per unit time that the island has
#' existed. This approach allows for comparing diversification
#' dynamics across islands of different geological ages.Events that predate
#' island formation are excluded from the count before computing rates.
#'
#' @param phy The matched phylogenetic tree.
#' @param events The data frame returned by \code{\link{age_correct_events}},
#'   which must include \code{node_age}, \code{island_age}, and
#'   \code{predates_island} columns.
#' @param island_ages A named numeric vector of island formation ages, in
#'   the same units as branch lengths.
#'
#' @return A data frame with one row per island and columns:
#'   \describe{
#'     \item{\code{island}}{Island name.}
#'     \item{\code{n_insitu}}{Number of valid in-situ events (excluding those
#'       predating the island).}
#'     \item{\code{island_age}}{Island formation age.}
#'     \item{\code{insitu_rate_corrected}}{In-situ events per unit time since
#'       island formation.}
#'   }
#'
#' @export
insitu_rate_by_island_age <- function(phy, events, island_ages) {
  if (!all(c("predates_island", "island_age") %in% names(events)))
    stop("Run age_correct_events() on your events table first.")

  # Exclude events predating island formation
  valid_events <- events[!events$predates_island & events$in_situ == TRUE, ]

  islands <- names(island_ages)

  out <- lapply(islands, function(isl) {
    n_insitu  <- sum(valid_events$island == isl, na.rm = TRUE)
    isl_age   <- island_ages[isl]
    data.frame(
      island = isl,
      n_insitu = n_insitu,
      island_age = isl_age,
      insitu_rate_corrected = if (isl_age > 0) n_insitu / isl_age else NA,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}
