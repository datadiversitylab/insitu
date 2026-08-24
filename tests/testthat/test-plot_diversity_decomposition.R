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

# Diversity decomposition
div <- decompose_diversity(phy = matched$phy,
                           events = map_events,
                           PAM = matched$PAM)
########
# This function's purpose is to plot, but we can test some elements to verify
#  accuracy without using snapshot or SVG tests

test_that("The correct values are plotted for Cuba", {
  p <- plot_diversity_decomposition(div)

  expect_equal(unname(p[,2]), c(10, 0, 1))
})

test_that("The correct values are plotted for Jamaica", {
  p <- plot_diversity_decomposition(div)

  expect_equal(unname(p[,5]), c(4, 1, 0))
})

