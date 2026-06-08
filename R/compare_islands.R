#' Statistical comparison of diversity decomposition metrics across multiple island systems
#'
#' Tests whether in-situ speciation, import, and export counts differ across
#' islands using a chi-square test, and whether in-situ and colonization rates
#' differ across islands using a Kruskal-Wallis test. Rates are computed
#' internally from the counts and total species per island.
#'
#' @param div A data frame returned by \code{decompose_diversity}.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{counts_test}}{Chi-square test comparing n_insitu, n_import,
#'       and n_export across islands.}
#'     \item{\code{insitu_rate_test}}{Kruskal-Wallis test on per-island
#'       in-situ speciation rates (n_insitu / n_species).}
#'     \item{\code{colonization_rate_test}}{Kruskal-Wallis test on per-island
#'       colonization rates ((n_import + n_export) / 2 / n_species).}
#'   }
#' @export
compare_islands <- function(div) {
  # counts comparison across islands
  counts_mat   <- as.matrix(div[, c("n_insitu", "n_import", "n_export")])
  rownames(counts_mat) <- div$island
  counts_test  <- chisq.test(counts_mat)

  # per-island rates
  insitu_rate       <- div$n_insitu / div$n_species
  colonization_rate <- ((div$n_import + div$n_export) / 2) / div$n_species

  insitu_rate_test       <- kruskal.test(insitu_rate       ~ div$island)
  colonization_rate_test <- kruskal.test(colonization_rate ~ div$island)

  list(
    counts_test            = counts_test,
    insitu_rate_test       = insitu_rate_test,
    colonization_rate_test = colonization_rate_test
  )
}
