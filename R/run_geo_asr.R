#' Runs geographic state ASR (discrete: island; mainland?)
#'
#' @param phy The phylogenetic tree associated with your data
#' @param PAM A presence-absence matrix reflecting where each species of
#' interest is located. The `match_island_phylo` function will create this PAM,
#' but if the user would prefer to input a custom PAM, make sure that it has a
#' column titled "locale" with the name of an island in each row, and each
#' subsequent column is titled with a species name. Each species column should
#' include either a 0 (absence) or a 1 (presence), signifying whether that
#' species occurs on the island in a given row.
#' @return A list of island-specific ancestral state reconstructions
#'
#' @export

run_geo_asr <- function(phy, PAM){
  # First, run an ancestral state reconstruction per island separately
  # Get the list of islands
  islands <- PAM$locale

  # Create an empty list to save reconstructions
  anc_recons <- list()

  # For each island in the PAM, run ace
  for(i in 1:length(islands)){
    # The full row except the first column (where locales are)
    # Needs to be a vector for ace
    island_states <- unlist(PAM[i, -1])
    # Save the reconstruction in the list under the name of the island
    anc_recons[[islands[i]]] <- ape::ace(x = island_states,
                                         phy = phy,
                                         type = "discrete",
                                         model = "ER")
  }

  # Return the list
  return(anc_recons)
}
