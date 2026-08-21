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

########

# Correct output class (data.frame)
test_that("map_insitu_events returns a data.frame", {
  map_events <- map_insitu_events(recons = asr_insitu,
                                  phy = matched$phy,
                                  PAM = matched$PAM,
                                  threshold = 0.5)
  expect_s3_class(map_events, "data.frame")
})

# Correct number of columns in the data.frame (6)
test_that("The resulting data.frame has the correct number of columns", {
  map_events <- map_insitu_events(recons = asr_insitu,
                                  phy = matched$phy,
                                  PAM = matched$PAM,
                                  threshold = 0.5)
  expect_equal(ncol(map_events), 6)
})
