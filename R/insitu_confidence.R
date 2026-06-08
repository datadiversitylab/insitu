#' Computes per-node support for in-situ event classification
#'
#' For each internal node, smmarizes the probability of in-situ
#' speciation across all islands by extracting the ancestor presence
#' probability directly from reconstructions.
#'
#' @param recons A list of ancestral state reconstructions returned by
#'   \code{run_geo_asr}.
#' @param events The data frame returned by \code{map_insitu_events}.
#'
#' @return A data frame with one row per node and columns:
#'   \describe{
#'     \item{\code{node}}{Internal node number.}
#'     \item{\code{island}}{Island associated with the node.}
#'     \item{\code{ancestor_presence_prob}}{Posterior probability of the
#'       ancestor being present on the island.}
#'     \item{\code{in_situ}}{Whether in-situ speciation was inferred.}
#'     \item{\code{support}}{Qualitative support level: \code{"high"}
#'       (prob >= 0.75), \code{"moderate"} (prob >= 0.5), or
#'       \code{"low"} (prob < 0.5).}
#'   }
#' @export
insitu_confidence <- function(recons, events) {
  in_situ_events <- events[events$in_situ == TRUE, ]

  in_situ_events$support <- cut(
    in_situ_events$ancestor_presence_prob,
    breaks = c(-Inf, 0.5, 0.75, Inf),
    labels = c("low", "moderate", "high"),
    right  = FALSE
  )

  in_situ_events[, c("node", "island", "ancestor_presence_prob",
                     "in_situ", "support")]
}
