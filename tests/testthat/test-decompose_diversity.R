# Read in example data and tree
# Read tree trimmed from Patton et al. 2021
tree <- ape::read.tree(system.file("extdata",
                                   "Patton_etal_trimmed.tree",
                                   package = "insitu"))

# Read species location dataframe
dat <- read.csv(system.file("extdata",
                            "anolis_dat.csv",
                            package = "insitu"))

# Ensure that species represented in both the tree and data are the same
matched <- match_island_phylo(phy = tree, locs = dat)

# Run ASR
asr_insitu <- run_geo_asr(phy = matched$phy, PAM = matched$PAM)

# Create event object
map_events <- map_insitu_events(recons = asr_insitu,
                                phy = matched$phy,
                                PAM = matched$PAM,
                                threshold = 0.5)

########

# Correct output class (data.frame)
test_that("decompose_diversity returns a data.frame", {
  div <- decompose_diversity(phy = matched$phy, events = map_events, PAM = matched$PAM)
  expect_s3_class(div, "data.frame")
})

# Correct number of columns (6)
test_that("decompose_diversity returns a data.frame with 6 columns", {
  div <- decompose_diversity(phy = matched$phy, events = map_events, PAM = matched$PAM)
  expect_equal(ncol(div), 6)
})

# Correct number of rows (equal to the number of islands)
test_that("decompose_diversity returns a data.frame with a row for each island", {
  div <- decompose_diversity(phy = matched$phy, events = map_events, PAM = matched$PAM)
  expect_equal(nrow(div), length(unique(matched$PAM$locale)))
})
