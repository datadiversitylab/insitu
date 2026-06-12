#' Uses stochastic character mapping to propagate ASR uncertainty through speciation classifications. Posterior distribution of in-situ event counts
#'
#' When multiple trees are provided, the function propagates ASR uncertainty
#' from the stochastic mapping itself, and phylogenetic uncertainty from
#' the distribution of trees.
#'
#' @param phy The phylogenetic tree associated with your data
#' @param PAM A presence-absence matrix reflecting where each species of
#' interest is located. The `match_island_phylo` function will create this PAM,
#' but if the user would prefer to input a custom PAM, make sure that it has a
#' column titled "locale" with the name of an island in each row, and each
#' subsequent column is titled with a species name. Each species column should
#' include either a 0 (absence) or a 1 (presence), signifying whether that
#' species occurs on the island in a given row.
#' @param model The transition model to use in `phytools::make.simmap`. See
#' documentation for `ape::ace` for more information. Default: ER
#' @param nsim The number of simulations to use in `phytools::make.simmap`
#' Default: 10
#' @return A list of island-specific ancestral state reconstructions
#'
#' @export

simmap_insitu <- function(phy, PAM, model = "ER", nsim = 10){
  # Accept either a single tree or a list of trees
  if (inherits(phy, "phylo")) phy <- list(phy)

  # First, run an ancestral state reconstruction per island separately
  # Get the list of islands
  islands <- PAM$locale

  # Run stochastic maps for each tree
  island_accumulator <- vector("list", length(islands))
  names(island_accumulator) <- islands

  for (tree in phy) {
    for (i in seq_along(islands)) {
      island_states <- unlist(PAM[i, -1])

      # When trees differ in tips
      island_states <- island_states[names(island_states) %in% tree$tip.label]

      ace_mat <- summary(
        phytools::make.simmap(
          tree = tree,
          x = island_states,
          model = model,
          nsim = nsim,
          Q = "mcmc",
          message = FALSE
        )
      )$ace

      if (is.null(island_accumulator[[islands[i]]])) {
        island_accumulator[[islands[i]]] <- ace_mat
      } else {
        # Average over shared nodes
        shared <- intersect(
          rownames(island_accumulator[[islands[i]]]),
          rownames(ace_mat)
        )
        island_accumulator[[islands[i]]][shared, ] <-
          (island_accumulator[[islands[i]]][shared, ] + ace_mat[shared, ]) / 2
      }
    }
  }

  return(island_accumulator)
}
