#' Simulates island assemblages under expected colonization, speciation, and
#' extinction parameters
#'
#' Generates a dated phylogeny and presence-absence matrix under a
#' birth-death-colonization process with known parameters. This function
#' returns both the simulated data and the true event table.
#'
#' In this function, a mainland tree of \code{n_mainland} species is generated
#' under a birth-death process. Geographic states are then
#' generated along the tree using a continuous-time Markov model where
#' \code{colonization_rate} governs transitions from mainland to island.
#' \code{back_rate} governs the reverse. Extinction is applied by pruning a
#' random fraction of island tips. Yet, we acknowledge that extinction is not
#' always random. The resulting tip states define the PAM, and nodes where
#' both parent and child reconstruct as island are the true in-situ events.
#'
#' @param n_tips Total number of tips in the simulated tree.
#' @param birth_rate Per-lineage speciation rate.
#' @param death_rate Per-lineage extinction rate. Set to \code{0} for a
#'   pure-birth tree.
#' @param colonization_rate Rate of transitions from mainland to island per
#'   unit branch length.
#' @param back_rate Rate of transitions from island back to mainland per unit
#'   branch length. Default: \code{0} (no back-colonization).
#' @param island_extinction_fraction Proportion of island tips to remove at
#'   random, simulating island extinction. Default: \code{0}.
#' @param n_islands Number of distinct islands to simulate. When greater than
#'   1, island states are distributed evenly across island-assigned tips.
#'   Default: \code{1}.
#' @param seed Optional random seed for reproducibility.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{phy}}{The simulated (and possibly pruned) phylogeny.}
#'     \item{\code{PAM}}{A presence-absence matrix in the format expected by
#'       \code{run_geo_asr} and \code{map_insitu_events}.}
#'     \item{\code{true_events}}{A data frame with the same structure as the
#'       output of \code{map_insitu_events}, containing the ground-truth
#'       in-situ classifications.}
#'   }
#'
#' @export
simulate_island <- function(n_tips,
                             birth_rate,
                             death_rate = 0,
                             colonization_rate,
                             back_rate = 0,
                             island_extinction_fraction = 0,
                             n_islands = 1,
                             seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Simulate a birth-death tree
  phy <- ape::rphylo(n = n_tips,
                     birth = birth_rate,
                     death = death_rate,
                     fossils = FALSE)

  n_species <- ape::Ntip(phy)

  # Build transition rate matrix: mainland (1) <-> island (2)
  # States: 1 = mainland, 2 = island
  Q <- matrix(c(-colonization_rate,  colonization_rate,
                 back_rate, -back_rate),
              nrow = 2, byrow = TRUE)

  # Simulate geographic states along the tree
  # All lineages start on the mainland (state 1)
  tip_states <- ape::rTraitDisc(
    phy,
    model = Q,
    k = 2,
    states = c("mainland", "island"),
    root.value = 1
  )

  # Apply island extinction: remove a fraction of island tips at random
  island_tips <- names(tip_states)[tip_states == "island"]
  if (island_extinction_fraction > 0 && length(island_tips) > 0) {
    n_extinct <- round(length(island_tips) * island_extinction_fraction)
    extinct_tips <- sample(island_tips, n_extinct)
    phy <- ape::drop.tip(phy, extinct_tips)
    tip_states <- tip_states[phy$tip.label]
  }

  n_species <- ape::Ntip(phy)

  # Assign island-state tips to islands (distribute evenly when n_islands > 1)
  island_tip_names <- names(tip_states)[tip_states == "island"]
  if (n_islands > 1 && length(island_tip_names) > 0) {
    assignment <- rep_len(paste0("island_", seq_len(n_islands)),
                          length(island_tip_names))
    assignment <- sample(assignment)  # randomize which tips go where
    tip_islands <- stats::setNames(assignment, island_tip_names)
  } else {
    tip_islands <- stats::setNames(
      rep("island_1", length(island_tip_names)),
      island_tip_names
    )
  }

  # Build PAM
  all_islands <- unique(c(tip_islands, "Mainland"))
  all_islands <- all_islands[all_islands != "Mainland"]
  locales <- c(all_islands, "Mainland")

  PAM <- data.frame(locale = locales, stringsAsFactors = FALSE)
  for (sp in phy$tip.label) {
    col <- rep(0L, length(locales))
    if (tip_states[sp] == "island") {
      isl_idx     <- which(locales == tip_islands[sp])
      col[isl_idx] <- 1L
    } else {
      mld_idx     <- which(locales == "Mainland")
      col[mld_idx] <- 1L
    }
    PAM[[sp]] <- col
  }

  # Reconstruct true in-situ events from the simulated states
  true_events <- data.frame()
  islands <- all_islands

  for (node in (n_species + 1):(n_species + phy$Nnode)) {
    children <- phy$edge[phy$edge[, 1] == node, 2]
    if (length(children) < 2) next

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

    desc1_ranges <- sapply(islands, function(isl)
      any(tip_islands[tips_desc1] == isl, na.rm = TRUE))
    desc2_ranges <- sapply(islands, function(isl)
      any(tip_islands[tips_desc2] == isl, na.rm = TRUE))
    common_islands <- islands[desc1_ranges & desc2_ranges]

    if (length(common_islands) == 1) {
      true_events <- rbind(true_events, data.frame(
        node = node,
        island = common_islands,
        in_situ = TRUE,
        stringsAsFactors = FALSE
      ))
    }
  }

  list(phy = phy, PAM = PAM, true_events = true_events)
}
