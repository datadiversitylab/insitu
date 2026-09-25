#' Simulates island assemblages under known colonization, speciation, and
#' extinction parameters
#'
#' Generates a dated phylogeny and presence-absence matrix under a
#' birth-death-colonization process with known parameters, returning both the
#' simulated data and the true event table for use in
#' \code{\link{insitu_power}}.
#'
#' The simulation models three distinct processes. Mainland lineages colonize
#' the island system at rate \code{mainland_colonization_rate}: the number of
#' independent colonization events is drawn from a Poisson distribution with
#' mean equal to that rate multiplied by the total mainland branch length. When
#' only one event is drawn, all island species are monophyletic. When multiple
#' events are drawn, the island assemblage is polyphyletic. Once on the island
#' system, lineages disperse between islands at rate \code{inter_island_rate}.
#' Speciation and extinction proceed at island-specific rates that may differ
#' from the mainland.
#'
#' True in-situ events are nodes where both descendant clades share exactly one
#' island AND descend from the same mainland colonization event. Nodes spanning
#' two independent colonization-derived clades are not in-situ at the system
#' level even if their descendants share an island.
#'
#' @param n_tips Total number of tips in the simulated tree.
#' @param diversification_rate Net mainland diversification rate (birth - death).
#' @param extinction_fraction Mainland relative extinction (death / birth),
#'   between 0 and 1. Default: \code{0}.
#' @param diversification_rate_island Island net diversification rate.
#'   Defaults to \code{diversification_rate} when not set.
#' @param extinction_fraction_island Island relative extinction.
#'   Defaults to \code{extinction_fraction} when not set.
#' @param mainland_colonization_rate Rate of mainland-to-island colonization
#'   events per unit of total mainland branch length.
#' @param inter_island_rate Rate of inter-island dispersal events per unit of
#'   total island branch length.
#' @param n_islands Number of distinct islands. Default: \code{1}.
#' @param n_mainland Number of mainland tips. Default: \code{1}.
#' @param seed Optional random seed. Default: \code{NULL}.
#'
#' @return A named list with elements \code{phy}, \code{PAM}, and
#'   \code{true_events}.
#'
#' @export
simulate_island <- function(n_tips,
                            diversification_rate,
                            extinction_fraction = 0,
                            diversification_rate_island = NULL,
                            extinction_fraction_island = NULL,
                            mainland_colonization_rate,
                            inter_island_rate,
                            n_islands = 1,
                            n_mainland = 1,
                            seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  if (extinction_fraction < 0 || extinction_fraction >= 1)
    stop("'extinction_fraction' must be between 0 and 1.")

  if (is.null(diversification_rate_island))
    diversification_rate_island <- diversification_rate
  if (is.null(extinction_fraction_island))
    extinction_fraction_island  <- extinction_fraction

  if (extinction_fraction_island < 0 || extinction_fraction_island >= 1)
    stop("'extinction_fraction_island' must be between 0 and 1.")

  # Derive birth and death rates
  birth_mland  <- diversification_rate / (1 - extinction_fraction)
  death_mland  <- birth_mland * extinction_fraction
  birth_island <- diversification_rate_island / (1 - extinction_fraction_island)
  death_island <- birth_island * extinction_fraction_island

  n_island_tips <- n_tips - n_mainland
  if (n_island_tips < 2L)
    stop("n_mainland leaves fewer than 2 island tips. ",
         "Reduce n_mainland or increase n_tips.")

  # Get all tip descendants of a node using the edge matrix only
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

  # Scale a tree to a target root-to-tip depth
  scale_to <- function(tree, depth) {
    max_d <- max(ape::node.depth.edgelength(tree))
    tree$edge.length <- tree$edge.length / max_d * depth
    tree
  }

  # Simulate mainland tree first to get branch length for Poisson draw
  phy_mland_tmp  <- ape::multi2di(
    ape::rphylo(n = n_mainland + 1L, birth = birth_mland,
                death = death_mland, fossils = FALSE)
  )
  total_bl_mland <- sum(phy_mland_tmp$edge.length)
  n_col_mainland <- max(1L, stats::rpois(
    1L, mainland_colonization_rate * total_bl_mland
  ))

  # Ensure at least 2 island tips per colonization event
  n_col_mainland <- min(n_col_mainland, floor(n_island_tips / 2L))
  n_col_mainland <- max(1L, n_col_mainland)

  # Distribute island tips across colonization events (min 2 per event)
  if (n_col_mainland == 1L) {
    island_tips_per_col <- n_island_tips
  } else {
    base                <- rep(2L, n_col_mainland)
    extra               <- n_island_tips - sum(base)
    add                 <- tabulate(
      sample(n_col_mainland, extra, replace = TRUE),
      nbins = n_col_mainland
    )
    island_tips_per_col <- base + add
  }

  # Build the actual mainland tree with one placeholder tip per colonization
  # event (fix: n_mainland + n_col_mainland tips, not n_mainland + 1)
  phy_mland <- ape::multi2di(
    ape::rphylo(n     = n_mainland + n_col_mainland,
                birth = birth_mland,
                death = death_mland,
                fossils = FALSE)
  )
  phy_mland$tip.label <- c(
    paste0("spM_", seq_len(n_mainland)),
    paste0("island_crown_", seq_len(n_col_mainland))
  )

  target_depth <- max(ape::node.depth.edgelength(phy_mland))

  # Each colonization group k gets tip labels "spIk_N", an unique prefix per
  # group allows robust crown identification after bind.tree
  for (k in seq_len(n_col_mainland)) {
    n_isl_k   <- island_tips_per_col[k]
    phy_isl_k <- ape::multi2di(
      ape::rphylo(n = n_isl_k, birth = birth_island,
                  death = death_island, fossils = FALSE)
    )
    phy_isl_k$tip.label <- paste0("spI", k, "_", seq_len(n_isl_k))
    phy_isl_k <- scale_to(phy_isl_k, target_depth * 0.5)

    ph_idx    <- which(phy_mland$tip.label == paste0("island_crown_", k))
    if (length(ph_idx) == 0L) next
    phy_mland <- ape::multi2di(
      ape::bind.tree(phy_mland, phy_isl_k, where = ph_idx)
    )
  }

  phy <- suppressMessages(
    phytools::force.ultrametric(phy_mland, method = "extend", message = FALSE)
  )

  n_species <- ape::Ntip(phy)
  n_nodes   <- ape::Nnode(phy)
  total_bl  <- sum(phy$edge.length)

  # Identify tip groups by prefix — robust to bind.tree renumbering
  island_tips_all   <- phy$tip.label[grepl("^spI",  phy$tip.label)]
  mainland_tips_all <- phy$tip.label[grepl("^spM_", phy$tip.label)]

  if (length(island_tips_all) == 0L) {
    warning("No island tips in final tree.")
    return(NULL)
  }

  # Find the crown node for each colonization group using group-specific prefixes
  island_crowns <- vapply(seq_len(n_col_mainland), function(k) {
    group_tips <- phy$tip.label[grepl(paste0("^spI", k, "_"), phy$tip.label)]
    if (length(group_tips) >= 2L) return(ape::getMRCA(phy, group_tips))
    # Single tip: return its parent internal node
    tip_idx    <- which(phy$tip.label == group_tips[1L])
    parent_row <- which(phy$edge[, 2] == tip_idx)
    if (length(parent_row) == 0L) return(NA_integer_)
    phy$edge[parent_row, 1L]
  }, integer(1L))

  island_crowns <- unique(island_crowns[!is.na(island_crowns)])

  if (length(island_crowns) == 0L) {
    warning("Could not identify island crown nodes.")
    return(NULL)
  }

  # Track which colonization event each island tip belongs to
  tip_col_id <- setNames(integer(length(island_tips_all)), island_tips_all)
  for (k in seq_len(n_col_mainland)) {
    group_tips <- phy$tip.label[grepl(paste0("^spI", k, "_"), phy$tip.label)]
    tip_col_id[group_tips] <- k
  }

  # Valid non-root internal nodes with >= 2 children
  all_nodes <- seq(n_species + 2L, n_species + n_nodes)
  all_nodes <- all_nodes[sapply(all_nodes, function(nd)
    sum(phy$edge[, 1] == nd) >= 2L)]

  # Inter-island dispersal within each island crown
  col_nodes      <- island_crowns
  node_to_island <- setNames(
    rep("island_1", length(island_crowns)),
    as.character(island_crowns)
  )

  for (k in seq_along(island_crowns)) {
    crown      <- island_crowns[k]
    crown_tips <- tips_of(phy, crown)

    island_internal <- all_nodes[sapply(all_nodes, function(nd) {
      nd != crown && all(tips_of(phy, nd) %in% crown_tips)
    })]

    n_inter     <- min(
      max(0L, stats::rpois(1L, inter_island_rate * total_bl)),
      length(island_internal),
      n_islands - 1L
    )
    inter_nodes <- if (n_inter > 0L && length(island_internal) > 0L)
      sample(island_internal, n_inter)
    else
      integer(0)

    all_col_k        <- c(crown, inter_nodes)
    isl_labels       <- paste0("island_",
                               ((seq_along(all_col_k) - 1L) %% n_islands) + 1L)
    node_to_island_k <- setNames(isl_labels, as.character(all_col_k))

    col_nodes                           <- c(col_nodes, inter_nodes)
    node_to_island                      <- c(node_to_island,
                                             node_to_island_k[as.character(inter_nodes)])
    node_to_island[as.character(crown)] <- isl_labels[1L]
  }

  # Assign each tip to its nearest colonization ancestor
  root_nd <- n_species + 1L

  col_ancestor <- function(tip_idx) {
    if (phy$tip.label[tip_idx] %in% mainland_tips_all) return(NA_integer_)
    current <- tip_idx
    repeat {
      parent_row <- which(phy$edge[, 2] == current)
      if (length(parent_row) == 0L) return(NA_integer_)
      parent <- phy$edge[parent_row, 1L]
      if (parent %in% col_nodes) return(parent)
      if (parent == root_nd)     return(NA_integer_)
      current <- parent
    }
  }

  tip_col_ancestor <- sapply(seq_len(n_species), col_ancestor)
  island_tip_idx   <- which(!is.na(tip_col_ancestor))

  if (length(island_tip_idx) == 0L) {
    warning("No island tips generated.")
    return(NULL)
  }

  tip_islands <- stats::setNames(
    node_to_island[as.character(tip_col_ancestor[island_tip_idx])],
    phy$tip.label[island_tip_idx]
  )
  tip_islands <- tip_islands[!is.na(tip_islands)]

  if (length(tip_islands) == 0L) {
    warning("No island tips could be assigned to islands.")
    return(NULL)
  }

  # Build PAM
  all_islands <- sort(unique(tip_islands))
  locales     <- c(all_islands, "Mainland")
  PAM         <- data.frame(locale = locales, stringsAsFactors = FALSE)

  for (sp in phy$tip.label) {
    col <- rep(0L, length(locales))
    col[which(locales == if (sp %in% names(tip_islands))
      tip_islands[sp] else "Mainland")] <- 1L
    PAM[[sp]] <- col
  }

  # True in-situ events
  # A node is in-situ when both descendant clades:
  #   (1) share exactly one island, AND
  #   (2) descend from the same mainland colonization event
  all_tip_regions <- setNames(rep("Mainland", n_species), phy$tip.label)
  all_tip_regions[names(tip_islands)] <- tip_islands

  true_events <- data.frame()

  for (node in seq(n_species + 1L, n_species + ape::Nnode(phy))) {
    children <- phy$edge[phy$edge[, 1] == node, 2]
    if (length(children) < 2L) next

    tips1 <- if (children[1L] <= n_species)
      phy$tip.label[children[1L]]
    else
      tips_of(phy, children[1L])

    tips2 <- if (children[2L] <= n_species)
      phy$tip.label[children[2L]]
    else
      tips_of(phy, children[2L])

    regions1 <- unique(all_tip_regions[tips1])
    regions2 <- unique(all_tip_regions[tips2])
    common   <- setdiff(intersect(regions1, regions2), "Mainland")

    if (length(common) != 1L) next

    # Both clades must share the same mainland colonization origin
    isl1 <- intersect(tips1, island_tips_all)
    isl2 <- intersect(tips2, island_tips_all)
    if (length(isl1) == 0L || length(isl2) == 0L) next
    if (length(intersect(unique(tip_col_id[isl1]),
                         unique(tip_col_id[isl2]))) == 0L) next

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

  list(phy = phy, PAM = PAM, true_events = true_events)
}
