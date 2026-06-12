#' Prepares a phylogeny for use in the insitu pipeline
#'
#' Handles non-ultrametricity and incomplete taxon sampling. The tree must
#' be fully bifurcating.
#'
#' @param phy A \code{phylo} object.
#' @param force_ultrametric Logical. If \code{TRUE} and the tree is not
#'   ultrametric, it is forced ultrametric using
#'   \code{phytools::force.ultrametric}. Default: \code{TRUE}.
#' @param sampling_fraction The proportion of extant species included in the
#'   phylogeny (between 0 and 1). Used to correct branch length-based rate
#'   estimates. Default: \code{1} (complete sampling assumed).
#'
#' @return The prepared \code{phylo} object, with \code{sampling_fraction}
#'   stored as an attribute.
#'
#' @export
prep_phylo <- function(phy, force_ultrametric = TRUE, sampling_fraction = 1) {
  if (!inherits(phy, "phylo"))
    stop("'phy' must be a phylo object.")
  if (sampling_fraction <= 0 || sampling_fraction > 1)
    stop("'sampling_fraction' must be between 0 and 1.")

  if (!ape::is.ultrametric(phy, tol = 1e-6)) {
    if (force_ultrametric) {
      message("Tree is not ultrametric. Forcing ultrametric using ",
              "phytools::force.ultrametric(method = 'extend'). ",
              "Internal node ages are preserved.")
      phy <- phytools::force.ultrametric(phy, method = "extend")
    } else {
      warning("Tree is not ultrametric. Some functions in insitu may fail. ",
              "Set force_ultrametric = TRUE to resolve this automatically.")
    }
  }

  if (sampling_fraction < 1)
    message("Sampling fraction set to ", sampling_fraction, ". ",
            "Rate estimates will be corrected for incomplete sampling.")

  attr(phy, "sampling_fraction") <- sampling_fraction
  phy
}
