plot_simulated_island <- function(sim) {
  phy     <- sim$phy
  islands <- sort(unique(sim$PAM$locale[sim$PAM$locale != "Mainland"]))
  colors  <- setNames(palette.colors(length(islands), "Okabe-Ito"), islands)

  # Assign tip colors
  tip_colors <- sapply(phy$tip.label, function(sp) {
    isl <- sim$PAM$locale[sim$PAM[[sp]] == 1]
    if (isl == "Mainland") "grey60" else colors[isl]
  })

  plot(phy,
       tip.color = tip_colors,
       main      = "Simulated island assemblage",
       cex       = 0.7)

  # Mark true in-situ nodes
  if (nrow(sim$true_events) > 0) {
    node_islands <- sim$true_events$island
    ape::nodelabels(node = sim$true_events$node,
                    pch  = 19,
                    col  = colors[node_islands],
                    cex  = 0.8)
  }

  legend("bottomleft",
         legend = c(islands, "Mainland", "in-situ node"),
         col    = c(colors, "grey60", "black"),
         pch    = c(rep(15, length(islands)), 15, 19),
         bty    = "n",
         cex    = 0.7)
}
