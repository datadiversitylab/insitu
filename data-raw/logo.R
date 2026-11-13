# Hex logo for insitu
# Run from the package root: source("data-raw/logo.R")

# Pointy-top hexagon
hex <- function(r = 1) {
  a <- seq(30, 330, by = 60) * pi / 180
  list(x = r * cos(a), y = r * sin(a))
}

# Filled circle, or a two-state pie when p is between 0 and 1
pie_at <- function(x, y, p, r, col_main, col_isl) {
  arc <- function(from, to) {
    t <- seq(from, to, length.out = 60) * pi / 180
    list(x = c(x, x + r * cos(t)), y = c(y, y + r * sin(t)))
  }
  if (p <= 0 || p >= 1) {
    polygon(arc(0, 360), col = if (p >= 1) col_isl else col_main, border = NA)
  } else {
    polygon(arc(90 - 360 * p, 90), col = col_isl, border = NA)
    polygon(arc(90, 90 + 360 * (1 - p)), col = col_main, border = NA)
  }
}

# Five-tip phylogeny drawn as segments (root at the bottom)
tree <- data.frame(
  x0 = c(-0.40, -0.20, -0.40, -0.30, 0.15, -0.30, -0.075, 0.15, 0.00, 0.00, 0.30, 0.20, 0.40, 0.20),
  y0 = c(0.22, 0.22, 0.22, -0.08, -0.08, -0.08, -0.20, 0.03, 0.14, 0.14, 0.14, 0.36, 0.36, 0.36),
  x1 = c(-0.40, -0.20, -0.20, -0.30, 0.15, 0.15, -0.075, 0.15, 0.00, 0.30, 0.30, 0.20, 0.40, 0.40),
  y1 = c(0.55, 0.55, 0.22, 0.22, 0.03, -0.08, -0.08, 0.14, 0.55, 0.14, 0.36, 0.55, 0.55, 0.36),
  island = rep(c(FALSE, TRUE), c(7, 7))
)

# Internal nodes with the probability of the island state
nodes <- data.frame(
  x = c(-0.075, -0.30, 0.15, 0.30),
  y = c(-0.08, 0.22, 0.14, 0.36),
  p = c(0.35, 0.05, 0.75, 0.95)
)

tips <- data.frame(x = c(-0.40, -0.20, 0.00, 0.20, 0.40), island = c(FALSE, FALSE, TRUE, TRUE, TRUE))

pal <- list(bg = "#1E3A4C", border = "#142833", main = "#E8D9B5", isl = "#3FB8A9", text = "#FFFFFF")

# Draws on the open device, whatever its format
draw_logo <- function() {
  par(mar = c(0, 0, 0, 0), lend = "square", ljoin = "mitre")
  plot.new()
  plot.window(xlim = c(-sqrt(3) / 2, sqrt(3) / 2), ylim = c(-1, 1), asp = 1, xaxs = "i", yaxs = "i")

  # Border is the outer hexagon showing around a slightly smaller inner one
  polygon(hex(1), col = pal$border, border = NA)
  polygon(hex(0.94), col = pal$bg, border = NA)

  # Mainland branches are drawn first so island branches sit on top where they meet
  segments(tree$x0, tree$y0, tree$x1, tree$y1, col = ifelse(tree$island, pal$isl, pal$main), lwd = 4)

  for (i in seq_len(nrow(nodes))) {
    pie_at(nodes$x[i], nodes$y[i], nodes$p[i], 0.058, pal$main, pal$isl)
  }
  for (i in seq_len(nrow(tips))) {
    pie_at(tips$x[i], 0.55, as.numeric(tips$island[i]), 0.04, pal$main, pal$isl)
  }

  text(0, -0.5, "insitu", col = pal$text, cex = 2.1, font = 2)
}

w <- 1.73
h <- 2
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)

png("man/figures/logo.png", width = w, height = h, units = "in", res = 300, bg = "transparent")
draw_logo()
dev.off()

# Cairo embeds the font in the PDF and is needed for SVG
if (!capabilities("cairo")) stop("PDF and SVG need cairo support in this R build")

cairo_pdf("man/figures/logo.pdf", width = w, height = h, bg = "transparent")
draw_logo()
dev.off()

svg("man/figures/logo.svg", width = w, height = h, bg = "transparent")
draw_logo()
dev.off()
