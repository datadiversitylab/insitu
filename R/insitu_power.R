#' Power analysis for in-situ speciation event recovery
#'
#' Runs repeated simulations under known colonization, speciation, and
#' extinction parameters and measures how reliably the insitu pipeline
#' recovers the true in-situ events. For each simulation, a tree and PAM
#' are generated using \code{\link{simulate_island}}, the full pipeline is
#' applied, and recovered events are compared to the known ground truth.
#'
#' Sensitivity measures the proportion of true in-situ events that the
#' pipeline correctly identifies. Precision measures the proportion of
#' pipeline in-situ calls that are truly in-situ. Both vary with tree size,
#' extinction rate, and the ASR threshold, so running this function across
#' a range of parameter values directly answers how trustworthy a given set
#' of in-situ classifications is likely to be.
#'
#' @param n_sim Number of simulation replicates. Default: \code{100}.
#' @param n_tips Total number of tips per simulated tree.
#' @param diversification_rate Net diversification rate (birth - death).
#' @param extinction_fraction Relative extinction rate (death / birth),
#'   between 0 and 1. Default: \code{0}.
#' @param diversification_rate_island Island-specific net diversification rate.
#'   Defaults to \code{diversification_rate} when not set.
#' @param extinction_fraction_island Island-specific relative extinction.
#'   Defaults to \code{extinction_fraction} when not set.
#' @param colonization_rate Rate of inter-island colonization events per unit
#'   of total branch length.
#' @param n_islands Number of islands to simulate. Default: \code{1}.
#' @param n_mainland Number of mainland tips. Default: \code{1}.
#' @param monophyletic_island Logical. Default: \code{TRUE}.
#' @param threshold ASR probability threshold passed to
#'   \code{map_insitu_events}. Default: \code{0.5}.
#' @param model Transition model passed to \code{run_geo_asr} or
#'   \code{simmap_insitu}. Default: \code{"ER"}.
#' @param use_simmap Logical. If \code{TRUE}, uses \code{simmap_insitu}
#'   instead of \code{run_geo_asr}. Default: \code{FALSE}.
#' @param nsim Number of stochastic maps per replicate when
#'   \code{use_simmap = TRUE}. Default: \code{10}.
#'
#' @return A data frame with one row per simulation replicate and columns
#'   \code{sim}, \code{n_true_insitu}, \code{n_recovered},
#'   \code{n_false_neg}, \code{n_false_pos}, \code{sensitivity}, and
#'   \code{precision}.
#'
#' @export
insitu_power <- function(n_sim                       = 100,
                          n_tips,
                          diversification_rate,
                          extinction_fraction         = 0,
                          diversification_rate_island = NULL,
                          extinction_fraction_island  = NULL,
                          colonization_rate,
                          n_islands                   = 1,
                          n_mainland                  = 1,
                          monophyletic_island         = TRUE,
                          threshold                   = 0.5,
                          model                       = "ER",
                          use_simmap                  = FALSE,
                          nsim                        = 10) {

  out <- lapply(seq_len(n_sim), function(s) {
    sim <- simulate_island(
      n_tips                      = n_tips,
      diversification_rate        = diversification_rate,
      extinction_fraction         = extinction_fraction,
      diversification_rate_island = diversification_rate_island,
      extinction_fraction_island  = extinction_fraction_island,
      colonization_rate           = colonization_rate,
      n_islands                   = n_islands,
      n_mainland                  = n_mainland,
      monophyletic_island         = monophyletic_island,
      seed                        = s
    )

    if (is.null(sim) || nrow(sim$true_events) == 0) {
      return(data.frame(sim           = s,
                        n_true_insitu = 0L,
                        n_recovered   = 0L,
                        n_false_neg   = 0L,
                        n_false_pos   = 0L,
                        sensitivity   = NA_real_,
                        precision     = NA_real_,
                        stringsAsFactors = FALSE))
    }

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
      map_insitu_events(recons    = recons,
                        phy       = sim$phy,
                        PAM       = sim$PAM,
                        threshold = threshold),
      error = function(e) NULL
    )
    if (is.null(events)) return(NULL)

    # Compare by node + island pair
    true_pairs     <- paste(sim$true_events$node,
                            sim$true_events$island, sep = "_")
    pipeline_pairs <- paste(events$node[events$in_situ == TRUE],
                            events$island[events$in_situ == TRUE], sep = "_")

    true_pos  <- length(intersect(pipeline_pairs, true_pairs))
    false_neg <- length(setdiff(true_pairs,       pipeline_pairs))
    false_pos <- length(setdiff(pipeline_pairs,   true_pairs))

    data.frame(
      sim           = s,
      n_true_insitu = length(true_pairs),
      n_recovered   = true_pos,
      n_false_neg   = false_neg,
      n_false_pos   = false_pos,
      sensitivity   = true_pos / max(1L, length(true_pairs)),
      precision     = true_pos / max(1L, length(pipeline_pairs)),
      stringsAsFactors = FALSE
    )
  })

  expected_cols <- c("sim", "n_true_insitu", "n_recovered",
                     "n_false_neg", "n_false_pos", "sensitivity", "precision")

  out_clean <- Filter(Negate(is.null), out)
  out_clean <- Filter(function(x) is.data.frame(x) &&
                        all(expected_cols %in% names(x)), out_clean)

  if (length(out_clean) == 0L) return(NULL)

  out_combined <- do.call(rbind, lapply(out_clean, function(x) x[expected_cols]))

  valid     <- out_combined[out_combined$n_true_insitu > 0, ]
  n_skipped <- nrow(out_combined) - nrow(valid)

  if (n_skipped > 0)
    message(n_skipped, " of ", n_sim, " replicates had no true in-situ events ",
            "and were excluded. Consider increasing colonization_rate or n_tips.")

  valid
}
