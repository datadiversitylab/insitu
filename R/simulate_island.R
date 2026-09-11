#' Simulates island assemblages under known colonization, speciation, and
#' extinction parameters
#'
#' Generates a dated phylogeny and presence-absence matrix under a
#' birth-death-colonization process with known parameters, returning both the
#' simulated data and the true event table for use in
#' \code{\link{insitu_power}}.
#'
#' When \code{monophyletic_island = TRUE}, all island species descend from a
#' single crown node, reflecting the biology of most island radiations where
#' all island species trace back to a single colonization from the mainland.
#' Inter-island colonization events are then placed within that clade. When
#' \code{monophyletic_island = FALSE}, colonization nodes are placed freely
#' across the tree.
#'
#' Speciation and extinction rates are specified as net diversification rate
#' and extinction fraction. Birth and death rates are derived internally as:
#' \code{birth = diversification / (1 - extinction_fraction)} and
#' \code{death = birth * extinction_fraction}.
#'
#' @param n_tips Total number of tips in the simulated tree.
#' @param diversification_rate Net diversification rate (birth - death).
#' @param extinction_fraction Relative extinction rate (death / birth),
#'   between 0 and 1. Default: \code{0}.
#' @param diversification_rate_island Island-specific net diversification rate.
#'   Defaults to \code{diversification_rate} when not set.
#' @param extinction_fraction_island Island-specific relative extinction.
#'   Defaults to \code{extinction_fraction} when not set.
#' @param colonization_rate Expected number of inter-island colonization events
#'   per unit of total branch length.
#' @param n_islands Number of distinct islands. Default: \code{1}.
#' @param n_mainland Number of tips forced to be mainland. Default: \code{1}.
#' @param monophyletic_island Logical. If \code{TRUE} (default), all island
#'   species descend from a single crown node. If \code{FALSE}, colonization
#'   nodes are placed freely across the tree.
#' @param seed Optional random seed. Default: \code{NULL}.
#'
#' @return A named list with elements \code{phy}, \code{PAM}, and
#'   \code{true_events}.
#'
#' @export
simulate_island <- function(n_tips,
                             diversification_rate,
                             extinction_fraction          = 0,
                             diversification_rate_island  = NULL,
                             extinction_fraction_island   = NULL,
                             colonization_rate,
                             n_islands                    = 1,
                             n_mainland                   = 1,
                             monophyletic_island          = TRUE,
                             seed                         = NULL) {

  if (!is.null(seed)) set.seed(seed)

  if (extinction_fraction < 0 || extinction_fraction >= 1)
    stop("'extinction_fraction' must be between 0 and 1.")

  # Default island rates to mainland rates if not specified
  if (is.null(diversification_rate_island))
    diversification_rate_island <- diversification_rate
  if (is.null(extinction_fraction_island))
    extinction_fraction_island <- extinction_fraction

  if (extinction_fraction_island < 0 || extinction_fraction_island >= 1)
    stop("'extinction_fraction_island' must be between 0 and 1.")

  # Derive birth and death rates
  birth_rate        <- diversification_rate  / (1 - extinction_fraction)
  death_rate        <- birth_rate * extinction_fraction
  birth_rate_island <- diversification_rate_island / (1 - extinction_fraction_island)
  death_rate_island <- birth_rate_island * extinction_fraction_island

  n_island_tips <- n_tips - n_mainland
  if (n_island_tips < 2L)
    stop("n_mainland leaves fewer than 2 island tips. ",
         "Reduce n_mainland or increase n_tips.")

  tips_of <- function(phy, nd) {
    n_sp  <- ape::Ntip(phy)
    desc  <- integer(0)
    queue <- nd
    while (length(queue) > 0) {
      children <- phy$edge[phy$edge[, 1] == queue[1], 2]
      desc     <- c(desc, children[children <= n_sp])
      queue    <- c(queue[-1], children[children > n_sp])
    }
    phy$tip.label[desc]
  }

  if (monophyletic_island) {
    phy_mland           <- ape::multi2di(
      ape::rphylo(n = n_mainland + 1L, birth = birth_rate,
                  death = death_rate, fossils = FALSE)
    )
    phy_mland$tip.label <- c(paste0("spM_", seq_len(n_mainland)), "island_clade")

    phy_island           <- ape::multi2di(
      ape::rphylo(n = n_island_tips, birth = birth_rate_island,
                  death = death_rate_island, fossils = FALSE)
    )
    phy_island$tip.label <- paste0("spI_", seq_len(n_island_tips))

    scale_to <- function(tree, depth) {
      max_d <- max(ape::node.depth.edgelength(tree))
      tree$edge.length <- tree$edge.length / max_d * depth
      tree
    }

    target_depth <- max(ape::node.depth.edgelength(phy_mland))
    phy_island   <- scale_to(phy_island, target_depth * 0.5)
    phy_mland    <- scale_to(phy_mland,  target_depth)

    placeholder_idx <- which(phy_mland$tip.label == "island_clade")
    phy <- ape::multi2di(
      ape::bind.tree(phy_mland, phy_island, where = placeholder_idx)
    )
    phy <- suppressMessages(phytools::force.ultrametric(phy, method = "extend", message = FALSE))

  } else {
    phy_mland           <- ape::multi2di(
      ape::rphylo(n = n_mainland, birth = birth_rate,
                  death = death_rate, fossils = FALSE)
    )
    phy_mland$tip.label  <- paste0("spM_", seq_len(n_mainland))

    phy_island            <- ape::multi2di(
      ape::rphylo(n = n_island_tips, birth = birth_rate_island,
                  death = death_rate_island, fossils = FALSE)
    )
    phy_island$tip.label  <- paste0("spI_", seq_len(n_island_tips))

    phy_mland$root.edge  <- 0.01
    phy_island$root.edge <- 0.01
    phy <- ape::multi2di(ape::bind.tree(phy_mland, phy_island))
    phy <- suppressMessages(phytools::force.ultrametric(phy, method = "extend"))
  }

  n_species <- ape::Ntip(phy)
  n_nodes   <- ape::Nnode(phy)
  total_bl  <- sum(phy$edge.length)

  island_tips_all   <- paste0("spI_", seq_len(n_island_tips))
  mainland_tips_all <- paste0("spM_", seq_len(n_mainland))

  all_nodes <- seq(n_species + 2L, n_species + n_nodes)
  all_nodes <- all_nodes[sapply(all_nodes, function(nd)
    sum(phy$edge[, 1] == nd) >= 2L)]

  island_crown <- ape::getMRCA(phy, island_tips_all)

  island_internal <- all_nodes[sapply(all_nodes, function(nd) {
    nd != island_crown &&
      all(tips_of(phy, nd) %in% island_tips_all)
  })]

  n_col     <- min(max(0L, stats::rpois(1, colonization_rate * total_bl)),
                   length(island_internal), n_islands - 1L)
  col_nodes <- if (n_col > 0L && length(island_internal) > 0L) {
    c(island_crown, sample(island_internal, n_col))
  } else {
    island_crown
  }

  node_to_island <- setNames(
    paste0("island_", ((seq_along(col_nodes) - 1L) %% n_islands) + 1L),
    as.character(col_nodes)
  )

  col_ancestor <- function(tip_idx) {
    if (phy$tip.label[tip_idx] %in% mainland_tips_all) return(NA_integer_)
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
  island_tip_idx   <- which(!is.na(tip_col_ancestor))

  if (length(island_tip_idx) == 0L) {
    warning("No island tips generated.")
    return(NULL)
  }


  tip_islands <- setNames(
    node_to_island[as.character(tip_col_ancestor[island_tip_idx])],
    phy$tip.label[island_tip_idx]
  )

  all_islands <- sort(unique(tip_islands))
  locales     <- c(all_islands, "Mainland")
  PAM         <- data.frame(locale = locales, stringsAsFactors = FALSE)

  for (sp in phy$tip.label) {
    col <- rep(0L, length(locales))
    col[which(locales == if (sp %in% names(tip_islands))
      tip_islands[sp] else "Mainland")] <- 1L
    PAM[[sp]] <- col
  }

  all_tip_regions[names(tip_islands)] <- tip_islands

  true_events <- data.frame()

  for (node in seq(n_species + 1L, n_species + ape::Nnode(phy))) {
    children <- phy$edge[phy$edge[, 1] == node, 2]
    if (length(children) < 2L) next

    regions1 <- if (children[1] <= n_species)
      all_tip_regions[phy$tip.label[children[1]]]
    else
      unique(all_tip_regions[tips_of(phy, children[1])])

    regions2 <- if (children[2] <= n_species)
      all_tip_regions[phy$tip.label[children[2]]]
    else
      unique(all_tip_regions[tips_of(phy, children[2])])

    common <- setdiff(intersect(regions1, regions2), "Mainland")

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
