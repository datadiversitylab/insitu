set.seed(42)

# Read in example tree
tree <- ape::read.tree(system.file("extdata/Patton_etal_trimmed.tree",
                                   package = "insitu"))

# Handle non-ultrametricity and incomplete taxon sampling
tree <- prep_phylo(tree)

# Read in example dataset
dat <- read.csv(system.file("extdata/anolis_dat.csv",
                            package = "insitu"))

# Ensure that species represented in both the tree and data are the same
matched <- match_island_phylo(phy = tree, locs = dat)

# Use the tree and PAM from the "matched" object to run an ancestral state reconstruction
asr_insitu <- run_geo_asr(phy = matched$phy, PAM = matched$PAM)

# Map events
map_events <- map_insitu_events(recons = asr_insitu,
                                phy = matched$phy,
                                PAM = matched$PAM,
                                threshold = 0.5)

########
# This function's purpose is to plot, but we can test some elements to verify
#  accuracy without using snapshot or SVG tests

test_that("The width of the plot is correct", {
  p <- plot_classified_phylo(phy = matched$phy,
                                events = map_events,
                                PAM = matched$PAM)
  expect_equal(round(p$rect$w, 2), 16.97)
})

test_that("The height of the plot is correct", {
  p <- plot_classified_phylo(phy = matched$phy,
                             events = map_events,
                             PAM = matched$PAM)
  expect_equal(round(p$rect$h, 2), 9.05)
})

test_that("The left margin of the plot is correct", {
  p <- plot_classified_phylo(phy = matched$phy,
                             events = map_events,
                             PAM = matched$PAM)
  expect_equal(round(p$rect$left, 2), 60.8)
})

test_that("The top margin of the plot is correct", {
  p <- plot_classified_phylo(phy = matched$phy,
                             events = map_events,
                             PAM = matched$PAM)
  expect_equal(round(p$rect$top, 2), 14.5)
})

