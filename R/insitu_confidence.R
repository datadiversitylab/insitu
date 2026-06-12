#' Computes per-node probabilities of island presence for all
#' internal nodes
#'
#' @param recons A list of ancestral state reconstructions returned by
#'   \code{run_geo_asr}.
#'
#' @return A data frame with one row per node per island and columns:
#'   \describe{
#'     \item{\code{node}}{Internal node number.}
#'     \item{\code{island}}{Island name.}
#'     \item{\code{prob}}{Probability of the ancestor being
#'       present on the island.}
#'     \item{\code{support}}{Qualitative support level: \code{"high"}
#'       (prob >= 0.75), \code{"moderate"} (prob >= 0.5), or
#'       \code{"low"} (prob < 0.5).}
#'   }
#' @export
insitu_confidence <- function(recons) {
  islands <- names(recons)
  results <- data.frame()

  for (island in islands) {
    lik <- if (typeof(recons[[1]]) == "list") {
      recons[[island]]$lik.anc
    } else {
      recons[[island]]
    }
    results <- rbind(results, data.frame(
      node = as.integer(rownames(lik)),
      island = island,
      prob = lik[, 2],
      stringsAsFactors = FALSE
    ))
  }

  results$support <- cut(
    results$prob,
    breaks = c(-Inf, 0.5, 0.75, Inf),
    labels = c("low", "moderate", "high"),
    right = FALSE
  )

  results
}
