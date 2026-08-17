#' Power analysis for in-situ speciation event recovery
#'
#' Runs repeated simulations under known colonization, speciation, and
#' extinction parameters and measures how reliably the insitu pipeline
#' recovers the true in-situ events. For each simulation, a tree and PAM
#' are generated using \code{\link{simulate_island}}, the full pipeline is
#' applied, and recovered events are compared to the known ground truth.
#'
#' Sensitivity measures the proportion of true in-situ events that the
#' pipeline correctly identifies. Specificity measures the proportion of
#' true non-in-situ nodes that are correctly left unclassified as in-situ.
#' Both vary with tree size, extinction rate, and the ASR threshold, so
#' running this function across a range of parameter values directly answers
#' how trustworthy a given set of in-situ classifications is likely to be.
#'
#' @param n_sim Number of simulation replicates. Default: \code{100}.
#' @param n_tips Total number of tips per simulated tree.
#' @param birth_rate Per-lineage speciation rate.
#' @param death_rate Per-lineage extinction rate. Default: \code{0}.
#' @param colonization_rate Rate of mainland-to-island transitions per unit
#'   branch length.
#' @param island_extinction_fraction Proportion of island tips to remove,
#'   simulating extinction. Default: \code{0}.
#' @param n_islands Number of islands to simulate. Default: \code{1}.
#' @param threshold ASR probability threshold passed to
#'   \code{map_insitu_events}. Default: \code{0.5}.
#' @param model Transition model passed to \code{run_geo_asr} or
#'   \code{simmap_insitu}. Default: \code{"ER"}.
#' @param use_simmap Logical. If \code{TRUE}, uses \code{simmap_insitu}
#'   instead of \code{run_geo_asr} to propagate ASR uncertainty through
#'   the power analysis. Default: \code{FALSE}.
#' @param nsim Number of stochastic maps per replicate when
#'   \code{use_simmap = TRUE}. Default: \code{10}.
#'
#' @return A data frame with one row per simulation replicate and columns:
#'   \describe{
#'     \item{\code{sim}}{Replicate number.}
#'     \item{\code{n_true_insitu}}{Number of true in-situ events in the
#'       simulated data.}
#'     \item{\code{n_recovered}}{Number of true in-situ events correctly
#'       recovered by the pipeline.}
#'     \item{\code{n_false_positive}}{Nodes classified as in-situ that were
#'       not true in-situ events.}
#'     \item{\code{sensitivity}}{True positive rate: \code{n_recovered /
#'       n_true_insitu}.}
#'     \item{\code{false_positive_rate}}{False positive rate: \code{
#'       n_false_positive / n_non_insitu_nodes}.}
#'   }
#'
#' @export
insitu_power <- function(n_sim = 100,
                          n_tips,
                          birth_rate,
                          death_rate = 0,
                          colonization_rate,
                          island_extinction_fraction = 0,
                          n_islands = 1,
                          threshold = 0.5,
                          model = "ER",
                          use_simmap = FALSE,
                          nsim = 10) {

  out <- lapply(seq_len(n_sim), function(s) {
    # Simulate
    sim <- simulate_island(
      n_tips = n_tips,
      birth_rate = birth_rate,
      death_rate = death_rate,
      colonization_rate = colonization_rate,
      island_extinction_fraction = island_extinction_fraction,
      n_islands = n_islands,
      seed = s
    )

    # Skip replicates where no island species were generated
    if (sum(sim$PAM[, -1]) == 0 || nrow(sim$true_events) == 0) {
      return(data.frame(sim = s,
                        n_true_insitu = 0,
                        n_recovered = 0,
                        n_false_positive = 0,
                        sensitivity = NA,
                        false_positive_rate = NA,
                        stringsAsFactors = FALSE))
    }

    # Run pipeline
    recons <- tryCatch(
      if (use_simmap) {
        simmap_insitu(phy = sim$phy, PAM = sim$PAM,
                      model = model, nsim = nsim)
      } else {
        run_geo_asr(phy = sim$phy, PAM = sim$PAM, model = model)
      },
      error = function(e) NULL
    )
    if (is.null(recons)) return(NULL)

    events <- tryCatch(
      map_insitu_events(recons = recons,
                        phy = sim$phy,
                        PAM = sim$PAM,
                        threshold = threshold),
      error = function(e) NULL
    )
    if (is.null(events)) return(NULL)

    # Compare recovered to true
    true_nodes <- sim$true_events$node
    recovered_nodes <- events$node[events$in_situ == TRUE]
    n_all_internal <- ape::Nnode(sim$phy)
    n_non_insitu <- n_all_internal - length(true_nodes)

    n_recovered <- length(intersect(true_nodes, recovered_nodes))
    n_false_positive <- length(setdiff(recovered_nodes, true_nodes))

    data.frame(
      sim = s,
      n_true_insitu = length(true_nodes),
      n_recovered = n_recovered,
      n_false_positive = n_false_positive,
      sensitivity = if (length(true_nodes) > 0)
        n_recovered / length(true_nodes) else NA,
      false_positive_rate = if (n_non_insitu > 0)
        n_false_positive / n_non_insitu else NA,
      stringsAsFactors    = FALSE
    )
  })

  do.call(rbind, Filter(Negate(is.null), out))
}
