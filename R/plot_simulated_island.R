#' Plots a simulated island assemblage with tips colored by island and
#' true in-situ nodes marked
#'
#' @param sim The list returned by \code{simulate_island}, containing
#'   \code{phy}, \code{PAM}, and \code{true_events}.
#'
#' @return Invisibly returns NULL.
#' @export
plot_simulated_island <- function(sim) {
  phy    <- sim$phy
  PAM    <- sim$PAM
  events <- sim$true_events

  islands <- PAM$locale[PAM$locale != "Mainland"]
  islands <- sort(unique(islands))
  colors  <- setNames(palette.colors(length(islands), "Okabe-Ito"), islands)

  # Assign a color to each tip
  tip_colors <- sapply(phy$tip.label, function(sp) {
    if (!sp %in% names(PAM)) return("grey60")
    isl <- PAM$locale[PAM[[sp]] == 1]
    if (length(isl) == 0 || isl == "Mainland") "grey60" else colors[isl]
  })

  plot(phy,
       tip.color = tip_colors,
       main      = "Simulated island assemblage",
       cex       = 0.7)

  if (nrow(events) > 0) {
    ape::nodelabels(node = events$node,
                    pch  = 19,
                    col  = colors[events$island],
                    cex  = 0.9)
  }

  legend("bottomleft",
         legend = c(islands, "Mainland", "in-situ node"),
         col    = c(colors, "grey60", "black"),
         pch    = c(rep(15, length(islands)), 15, 19),
         bty    = "n",
         cex    = 0.7)

  invisible(NULL)
}
