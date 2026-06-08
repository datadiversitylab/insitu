#' Flags internal nodes as within-island or between-island based on
#' descendant tip composition
#'
#' Classifies each internal node as within-island if all descendant tips
#' occur on the same island, between-island if descendants span multiple
#' islands, or ambiguous if no island assignment can be made. This
#' classification is based on tip composition with no ancestral
#' state reconstruction.
#'
#' @param phy The matched phylogenetic tree.
#' @param PAM A presence-absence matrix from \code{match_island_phylo}.
#'
#' @return A data frame with one row per internal node and columns:
#'   \describe{
#'     \item{\code{node}}{Internal node number.}
#'     \item{\code{status}}{One of \code{"within_island"},
#'       \code{"between_island"}, or \code{"ambiguous"}.}
#'     \item{\code{island}}{Island name for within-island nodes, \code{NA}
#'       otherwise.}
#'   }
#' @export
identify_island_nodes <- function(phy, PAM) {
  n_tips   <- ape::Ntip(phy)
  islands  <- PAM$locale
  results  <- data.frame()

  for (node in (n_tips + 1):(n_tips + phy$Nnode)) {
    # Collect all descendant tips
    desc  <- integer(0)
    queue <- node
    while (length(queue) > 0) {
      children <- phy$edge[phy$edge[, 1] == queue[1], 2]
      tips     <- children[children <= n_tips]
      internal <- children[children >  n_tips]
      desc     <- c(desc, tips)
      queue    <- c(queue[-1], internal)
    }

    tip_names <- phy$tip.label[desc]
    tip_names <- tip_names[tip_names %in% names(PAM)]

    if (length(tip_names) == 0L) {
      results <- rbind(results, data.frame(
        node   = node,
        status = "ambiguous",
        island = NA,
        stringsAsFactors = FALSE
      ))
      next
    }

    # Islands represented by descendant tips
    tip_islands <- lapply(tip_names, function(sp) islands[PAM[[sp]] == 1])
    all_islands <- unique(unlist(tip_islands))

    if (length(all_islands) == 1) {
      results <- rbind(results, data.frame(
        node   = node,
        status = "within_island",
        island = all_islands,
        stringsAsFactors = FALSE
      ))
    } else if (length(all_islands) > 1) {
      results <- rbind(results, data.frame(
        node   = node,
        status = "between_island",
        island = NA_character_,
        stringsAsFactors = FALSE
      ))
    } else {
      results <- rbind(results, data.frame(
        node   = node,
        status = "ambiguous",
        island = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
  }

  results
}
