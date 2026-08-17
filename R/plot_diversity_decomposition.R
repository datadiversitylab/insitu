#' Stacked bar chart of in-situ, import, and export diversity per island
#'
#' For each island, plots a stacked bar showing the number of in-situ
#' speciation, import, and export events as returned by
#' \code{decompose_diversity}.
#'
#' @param div A data frame returned by \code{decompose_diversity}.
#'
#' @return Invisibly returns the matrix used for plotting.
#' @export
plot_diversity_decomposition <- function(div) {
  mat <- rbind(
    in_situ = div$n_insitu,
    import  = div$n_import,
    export  = div$n_export
  )
  colnames(mat) <- div$island

  colors <- c(in_situ = "#E69F00", import = "#56B4E9", export = "#009E73")

  graphics::barplot(mat,
          col = colors,
          border  = NA,
          legend = rownames(mat),
          args.legend = list(bty = "n", cex = 0.8),
          ylab = "Number of events",
          las = 2,
          main = "Diversity decomposition per island")

  invisible(mat)
}
