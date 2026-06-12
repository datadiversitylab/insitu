#' Applies decision rules to reconstructed nodes to classify each speciation event as in-situ, colonization, or ambiguous
#'
#' @param recons A list of ancestral state reconstructions returned by `run_geo_asr`.
#' We do not recommend inputting a custom list for this function, but if you
#' need to, make sure that the list is composed of `ape::ace` objects, each
#' named after a locality. Refer to `run_geo_asr` for more information about
#' how the original function estimates ancestral states.
#' @param phy The phylogenetic tree associated with your data
#' @param PAM A presence-absence matrix reflecting where each species of
#' interest is located. The `match_island_phylo` function will create this PAM,
#' but if the user would prefer to input a custom PAM, make sure that it has a
#' column titled "locale" with the name of an island in each row, and each
#' subsequent column is titled with a species name. Each species column should
#' include either a 0 (absence) or a 1 (presence), signifying whether that
#' species occurs on the island in a given row.
#' @param threshold The threshold of probability at which it is reasonable to
#' infer that a given ancestral node occurred on the same island as a given
#' set of tips. Default: 0.5
#' @return A dataframe including node numbers, islands associated with those
#' nodes, and whether an in situ speciation event was likely at that node/island
#' combination.
#'
#' @export

map_insitu_events <- function(recons, phy, PAM, threshold = 0.5){
  # Initialize result dataframe
  results <- data.frame()

  # Get the number of species
  n_species <- length(phy$tip.label)

  # Get the list of islands from the PAM
  islands <- PAM$locale

  # For each relevant node, look at descendants and determine whether they are
  #   on the same island. If so, refer to the ancestral state reconstruction to
  #   infer whether an in situ speciation event occurred on this island
  for(node in (n_species + 1):(n_species + phy$Nnode)){
    # Get direct descendant nodes
    children <- phy$edge[phy$edge[,1] == node, 2]

    # Get descendant tips per child node
    tips_desc1 <- if(children[1] <= length(phy$tip.label)){
      phy$tip.label[children[1]]
    } else {
      ape::extract.clade(phy, children[1])$tip.label
    }

    tips_desc2 <- if(children[2] <= length(phy$tip.label)){
      phy$tip.label[children[2]]
    } else {
      ape::extract.clade(phy, children[2])$tip.label
    }

    # Species distributions for descendants
    # (TRUE/FALSE in each island from PAM rows)
    desc1_ranges <- rowSums(PAM[,tips_desc1, drop = F]) > 0
    desc2_ranges <- rowSums(PAM[,tips_desc2, drop = F]) > 0

    # Identify common islands between descendants
    common_islands <- islands[desc1_ranges & desc2_ranges]

    # Infer in situ speciation event only if exactly one common island
    if(length(common_islands) == 1){
      island <- common_islands

      # IF USED run_geo_asr
      if(typeof(recons[[1]]) == "list"){
        # Probability ancestor present on the same island
        anc_prob <- recons[[island]]$lik.anc[node - n_species, 2]

        # Confirm ancestor present with high likelihood
        in_situ <- anc_prob >= threshold

        results <- rbind(results, data.frame(node = node,
                                             island = island,
                                             ancestor_presence_prob = anc_prob,
                                             in_situ = in_situ))
      } else {
        # IF USED simmap_insitu
        # Probability ancestor present on the same island
        anc_prob <- recons[[island]][node - n_species, 2]

        # Confirm ancestor present with high likelihood
        in_situ <- anc_prob >= 0.5

        results <- rbind(results, data.frame(node = node,
                                             island = island,
                                             ancestor_presence_prob = anc_prob,
                                             in_situ = in_situ))
      }
    }
    # else ambiguous (multiple or zero islands), do not infer
  } # End for

  # Detect transitions: best island per node above threshold
  best_region <- function(nd) {
    if (nd <= n_species) {
      tip_name    <- phy$tip.label[nd]
      if (!tip_name %in% names(PAM)) return(NA)
      tip_islands <- islands[PAM[[tip_name]] == 1]
      if (length(tip_islands) == 1L) tip_islands else NA
    } else {
      probs <- sapply(islands, function(isl) {
        if (typeof(recons[[1]]) == "list") recons[[isl]]$lik.anc[nd - n_species, 2]
        else recons[[isl]][nd - n_species, 2]
      })
      if (max(probs) >= threshold) islands[which.max(probs)] else NA
    }
  }

  transitions <- data.frame()

  # For internal nodes
  for (i in seq_len(nrow(phy$edge))) {
    parent_nd <- phy$edge[i, 1]
    child_nd <- phy$edge[i, 2]
    parent_region <- best_region(parent_nd)
    child_region  <- best_region(child_nd)
    if (is.na(parent_region) || is.na(child_region)) next
    if (parent_region == child_region) next
    transitions <- rbind(transitions, data.frame(
      node = child_nd,
      export = parent_region,
      import = child_region,
      stringsAsFactors = FALSE
    ))
  }

  # For tips
  for (tip_idx in seq_len(n_species)) {
    parent_nd     <- phy$edge[phy$edge[, 2] == tip_idx, 1]
    parent_region <- best_region(parent_nd)
    tip_name      <- phy$tip.label[tip_idx]
    if (!tip_name %in% names(PAM)) next
    tip_islands   <- islands[PAM[[tip_name]] == 1]
    for (tip_island in tip_islands) {
      if (is.na(parent_region)) next
      if (parent_region == tip_island) next
      transitions <- rbind(transitions, data.frame(
        node = tip_idx,
        export = parent_region,
        import = tip_island,
        stringsAsFactors = FALSE
      ))
    }
  }

  results$export <- NA
  results$import <- NA
  transitions$ancestor_presence_prob <- NA
  transitions$in_situ <- NA
  transitions$island <- "TRANSITION"
  results <- rbind(results, transitions)
  return(results)
}
