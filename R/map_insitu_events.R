#' Applies decision rules to reconstructed nodes to classify each speciation event as in-situ, colonization, or ambiguous
#'
#' @param recons A list of ancestral state reconstructions returned by \code{run_geo_asr}.
#'   We do not recommend inputting a custom list for this function, but if you
#'   need to, make sure that the list is composed of \code{ape::ace} objects, each
#'   named after a locality. Refer to \code{run_geo_asr} for more information about
#'   how the original function estimates ancestral states.
#' @param phy The phylogenetic tree associated with your data.
#' @param PAM A presence-absence matrix reflecting where each species of
#'   interest is located. The \code{match_island_phylo} function will create this PAM,
#'   but if the user would prefer to input a custom PAM, make sure that it has a
#'   column titled "locale" with the name of an island in each row, and each
#'   subsequent column is titled with a species name. Each species column should
#'   include either a 0 (absence) or a 1 (presence), signifying whether that
#'   species occurs on the island in a given row.
#' @param threshold The threshold of probability at which it is reasonable to
#'   infer that a given ancestral node occurred on the same island as a given
#'   set of tips. Default: 0.5.
#'
#' @return A data frame including node numbers, islands associated with those
#'   nodes, whether an in-situ speciation event was likely at that node/island
#'   combination, and for transitions, the source island (export) and
#'   destination island (import).
#'
#' @export
map_insitu_events <- function(recons, phy, PAM, threshold = 0.5){
  results   <- data.frame()
  n_species <- length(phy$tip.label)
  islands   <- PAM$locale

  # helper: most likely island for an internal node above threshold
  parent_island_of <- function(nd) {
    if(length(nd) == 0 || nd <= n_species) return(NA)
    probs <- sapply(islands, function(isl) {
      if(typeof(recons[[1]]) == "list"){
        recons[[isl]]$lik.anc[nd - n_species, 2]
      } else {
        recons[[isl]][nd - n_species, 2]
      }
    })
    if(max(probs) >= threshold) islands[which.max(probs)] else NA
  }

  # internal nodes
  for(node in (n_species + 1):(n_species + phy$Nnode)){
    children <- phy$edge[phy$edge[, 1] == node, 2]

    tips_desc1 <- if(children[1] <= n_species){
      phy$tip.label[children[1]]
    } else {
      ape::extract.clade(phy, children[1])$tip.label
    }

    tips_desc2 <- if(children[2] <= n_species){
      phy$tip.label[children[2]]
    } else {
      ape::extract.clade(phy, children[2])$tip.label
    }

    desc1_ranges   <- rowSums(PAM[, tips_desc1, drop = FALSE]) > 0
    desc2_ranges   <- rowSums(PAM[, tips_desc2, drop = FALSE]) > 0
    common_islands <- islands[desc1_ranges & desc2_ranges]

    if(length(common_islands) == 1){
      island <- common_islands

      anc_prob <- if(typeof(recons[[1]]) == "list"){
        recons[[island]]$lik.anc[node - n_species, 2]
      } else {
        recons[[island]][node - n_species, 2]
      }

      in_situ       <- anc_prob >= threshold
      parent_nd     <- phy$edge[phy$edge[, 2] == node, 1]
      parent_island <- parent_island_of(parent_nd)
      export_island <- if(!is.na(parent_island) && parent_island != island) parent_island else NA
      import_island <- if(!is.na(parent_island) && parent_island != island) island else NA

      results <- rbind(results, data.frame(
        node                   = node,
        island                 = island,
        ancestor_presence_prob = anc_prob,
        in_situ                = in_situ,
        export                 = export_island,
        import                 = import_island,
        stringsAsFactors       = FALSE
      ))
    }
  }

  # tips: in-situ not possible, only import and export
  for(tip_idx in seq_len(n_species)){
    tip_name  <- phy$tip.label[tip_idx]
    if(!tip_name %in% names(PAM)) next
    tip_islands <- islands[PAM[[tip_name]] == 1]
    if(length(tip_islands) == 0L) next

    parent_nd     <- phy$edge[phy$edge[, 2] == tip_idx, 1]
    parent_island <- parent_island_of(parent_nd)
    if(is.na(parent_island)) next

    for(island in tip_islands){
      export_island <- if(parent_island != island) parent_island else NA
      import_island <- if(parent_island != island) island else NA

      results <- rbind(results, data.frame(
        node                   = tip_idx,
        island                 = island,
        ancestor_presence_prob = NA_real_,
        in_situ                = FALSE,
        export                 = export_island,
        import                 = import_island,
        stringsAsFactors       = FALSE
      ))
    }
  }

  results
}
