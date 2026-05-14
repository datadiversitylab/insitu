#' Uses stochastic character mapping to propagate ASR uncertainty through speciation classifications. Posterior distribution of in-situ event counts
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
  # First, run an ancestral state reconstruction per island separately
  # Get the list of islands
  islands <- PAM$locale

  # Create an empty list to save reconstructions
  sim_recons <- list()

  # For each island in the PAM, run make.simmap and save the output for
  #  each internal node
  for(i in 1:length(islands)){
    # The full row except the first column (where locales are)
    # Needs to be a vector for ace
    island_states <- unlist(PAM[i, -1])
    # Save the per-node reconstruction in the list under the name of the island
    sim_recons[[islands[i]]] <- summary(phytools::make.simmap(tree = phy,
                                                              x = island_states,
                                                              model = model,
                                                              nsim = nsim,
                                                              Q = "mcmc",
                                                              message = FALSE))$ace
  }

  # Return the list
  return(sim_recons)
}
