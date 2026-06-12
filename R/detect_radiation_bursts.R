#' Identifies nodes where multiple sequential in-situ speciation events occur
#' within a defined time window
#'
#' For each internal node marked as in-situ, counts how many additional in-situ
#' speciation events occur among its descendants within a given time window.
#' Each event represents a speciation that kept the lineage on the same island.
#' Clades exceeding the minimum event count are flagged as candidate adaptive
#' radiations.
#'
#' @param phy The matched phylogenetic tree.
#' @param events The data frame returned by \code{map_insitu_events}.
#' @param time_window The maximum time (Ma) span within which in-situ events must
#'   occur to be considered a burst.
#' @param min_events The minimum number of in-situ speciation events within
#'   the time window to qualify as a burst.
#' @param plot Logical. Whether to plot the phylogeny with burst nodes
#'   highlighted. Default: \code{TRUE}.
#'
#' @return A data frame with one row per detected burst and columns:
#'   \describe{
#'     \item{\code{node}}{Root node of the burst clade (MRCA).}
#'     \item{\code{island}}{Island on which the burst occurred.}
#'     \item{\code{n_insitu_events}}{Number of in-situ speciation events within
#'       the time window.}
#'     \item{\code{n_species}}{Total extant species descending from the burst clade.}
#'     \item{\code{burst_timeframe}}{Time elapsed from the MRCA to the youngest
#'       in-situ event in the burst.}
#'     \item{\code{time_window}}{The time window used.}
#'   }
#' @export
detect_radiation_bursts <- function(phy, events, time_window = 5,
                                    min_events = 3, plot = TRUE) {
  insitu_nodes <- events$node[events$in_situ == TRUE]
  insitu_nodes <- na.omit(insitu_nodes)
  node_ages <- ape::node.depth.edgelength(phy)
  n_tips <- ape::Ntip(phy)
  results <- data.frame()

  for (nd in insitu_nodes) {
    island <- events$island[events$node == nd & events$in_situ == TRUE][1]
    nd_age <- node_ages[nd]

    # Collect all descendants
    all_desc <- integer(0)
    queue    <- nd
    while (length(queue) > 0) {
      children <- phy$edge[phy$edge[, 1] == queue[1], 2]
      all_desc <- c(all_desc, children)
      queue <- c(queue[-1], children[children > n_tips])
    }

    # in-situ events among descendant internal nodes within time window
    desc_internal <- all_desc[all_desc > n_tips]
    desc_insitu <- intersect(desc_internal, insitu_nodes)
    desc_ages <- node_ages[desc_insitu]
    in_window <- desc_insitu[abs(desc_ages - nd_age) <= time_window]
    n_events <- 1L + length(in_window)

    if (n_events >= min_events) {
      n_species <- sum(all_desc <= n_tips)
      in_window_ages <- node_ages[in_window]
      burst_timeframe <- if (length(in_window) > 0)
        max(in_window_ages) - nd_age else 0

      results <- rbind(results, data.frame(
        node = nd,
        island = island,
        n_insitu_events = n_events,
        n_species = n_species,
        burst_timeframe = burst_timeframe,
        time_window = time_window,
        stringsAsFactors = FALSE
      ))
    }
  }

  if (plot && nrow(results) > 0) {
    islands <- unique(results$island)
    colors  <- setNames(palette.colors(length(islands), "Okabe-Ito"), islands)

    plot(phy,
         main = expression(paste("Adaptive radiation bursts")),
         cex  = 0.8)

    burst_nodes <- results$node
    burst_islands <- results$island

    ape::nodelabels(node = burst_nodes, pch = 19,
                    col  = colors[burst_islands], cex = 1.5)
    ape::nodelabels(text  = burst_islands,
                    node  = burst_nodes,
                    cex   = 0.8,
                    frame = "none",
                    adj   = c(-0.2, -0.4))

    legend("bottomleft",
           legend = islands,
           pch    = 19,
           col    = colors,
           bty    = "n",
           cex    = 0.7)
  }

  return(results)
}
