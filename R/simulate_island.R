#' Simulates island assemblages under known colonization, speciation, and
#' extinction parameters
#'
#' Generates a dated phylogeny and presence-absence matrix under a
#' birth-death-colonization process with known parameters, returning both the
#' simulated data and the true event table for use in
#' \code{\link{insitu_power}}.
#'
#' The simulation places colonization events directly on the tree by randomly
#' selecting internal nodes as colonization origins. The number of colonization
#' events is drawn from a Poisson distribution with mean equal to
#' \code{colonization_rate} multiplied by the total branch length of the tree.
#' Tips descending from a colonization node are assigned to an island; all
#' other tips are mainland. True in-situ events are nodes where both descendant
#' clades share exactly one island, computed directly from the simulated
#' topology without any ASR.
#'
#' @param n_tips Total number of tips in the simulated tree.
#' @param birth_rate Per-lineage speciation rate.
#' @param death_rate Per-lineage extinction rate. Default: \code{0}.
#' @param colonization_rate Expected number of colonization events per unit
#'   of total branch length. Controls how many internal nodes are selected
#'   as colonization origins.
#' @param island_extinction_fraction Proportion of island tips to remove at
#'   random, simulating extinction. Default: \code{0}.
#' @param n_islands Number of distinct islands. Colonization nodes are
#'   distributed across islands in sequence. Default: \code{1}.
#' @param seed Optional random seed. Default: \code{NULL}.
#'
#' @return A named list with elements \code{phy}, \code{PAM}, and
#'   \code{true_events}. \code{true_events} has the same structure as the
#'   output of \code{map_insitu_events} and serves as the ground truth for
#'   power analysis.
#'
#' @export
simulate_island <- function(n_tips,
                            birth_rate,
                            death_rate                 = 0,
                            colonization_rate,
                            island_extinction_fraction = 0,
                            n_islands                  = 1,
                            seed                       = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Simulate birth-death tree
  phy          <- ape::rphylo(n = n_tips, birth = birth_rate,
                              death = death_rate, fossils = FALSE)
  phy$tip.label <- paste0("sp_", seq_len(ape::Ntip(phy)))
  n_species    <- ape::Ntip(phy)
  n_nodes      <- ape::Nnode(phy)

  # Draw number of colonization events from Poisson
  # Exclude root so mainland outgroup is always present
  total_bl     <- sum(phy$edge.length)
  n_col        <- max(1L, stats::rpois(1, lambda = colonization_rate * total_bl))
  candidates   <- seq(n_species + 2L, n_species + n_nodes)

  if (length(candidates) == 0) {
    warning("Tree has no internal nodes available for colonization. ",
            "Try increasing n_tips.")
    return(NULL)
  }

  n_col     <- min(n_col, length(candidates))
  col_nodes <- sample(candidates, n_col)

  # For each tip, find its closest colonization ancestor (if any)
  col_ancestor <- function(tip_idx) {
    current <- tip_idx
    repeat {
      parent_row <- which(phy$edge[, 2] == current)
      if (length(parent_row) == 0L) return(NA_integer_)
      parent <- phy$edge[parent_row, 1L]
      if (parent %in% col_nodes) return(parent)
      current <- parent
    }
  }

  tip_col_ancestor <- sapply(seq_len(n_species), col_ancestor)

  island_tip_idx <- which(!is.na(tip_col_ancestor))
  if (length(island_tip_idx) == 0L) {
    warning("No island tips generated.")
    return(NULL)
  }

  # Assign colonization nodes to islands (cycle through islands)
  node_to_island <- stats::setNames(
    paste0("island_", ((seq_along(col_nodes) - 1L) %% n_islands) + 1L),
    as.character(col_nodes)
  )

  tip_islands <- stats::setNames(
    node_to_island[as.character(tip_col_ancestor[island_tip_idx])],
    phy$tip.label[island_tip_idx]
  )

  # Apply island extinction
  if (island_extinction_fraction > 0 && length(tip_islands) > 0) {
    n_extinct <- round(length(tip_islands) * island_extinction_fraction)
    if (n_extinct > 0L) {
      extinct      <- sample(names(tip_islands), n_extinct)
      phy          <- ape::drop.tip(phy, extinct)
      tip_islands  <- tip_islands[!names(tip_islands) %in% extinct]
      n_species    <- ape::Ntip(phy)
    }
  }

  # Build PAM
  all_islands <- sort(unique(tip_islands))
  locales     <- c(all_islands, "Mainland")
  PAM         <- data.frame(locale = locales, stringsAsFactors = FALSE)

  for (sp in phy$tip.label) {
    col <- rep(0L, length(locales))
    if (sp %in% names(tip_islands)) {
      col[which(locales == tip_islands[sp])] <- 1L
    } else {
      col[which(locales == "Mainland")] <- 1L
    }
    PAM[[sp]] <- col
  }

  # Region lookup for all tips
  all_tip_regions <- stats::setNames(rep("Mainland", n_species), phy$tip.label)
  all_tip_regions[names(tip_islands)] <- tip_islands

  # Compute true in-situ events directly from topology
  true_events <- data.frame()

  for (node in seq(n_species + 1L, n_species + ape::Nnode(phy))) {
    children <- phy$edge[phy$edge[, 1] == node, 2]
    if (length(children) < 2L) next

    tips_desc1 <- if (children[1] <= n_species) {
      phy$tip.label[children[1]]
    } else {
      ape::extract.clade(phy, children[1])$tip.label
    }

    tips_desc2 <- if (children[2] <= n_species) {
      phy$tip.label[children[2]]
    } else {
      ape::extract.clade(phy, children[2])$tip.label
    }

    common <- intersect(
      unique(all_tip_regions[tips_desc1]),
      unique(all_tip_regions[tips_desc2])
    )
    common <- common[common != "Mainland"]

    if (length(common) == 1L) {
      true_events <- rbind(true_events, data.frame(
        node                   = node,
        island                 = common,
        ancestor_presence_prob = NA_real_,
        in_situ                = TRUE,
        export                 = NA_character_,
        import                 = NA_character_,
        stringsAsFactors       = FALSE
      ))
    }
  }

  list(phy = phy, PAM = PAM, true_events = true_events)
}
