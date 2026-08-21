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
test_that("insitu_confidence returns a data.frame", {
  conf <- insitu_confidence(asr_insitu)
  expect_s3_class(conf, "data.frame")
})

# Correct number of columns (4)
test_that("insitu_confidence returns a data.frame with 6 columns", {
  conf <- insitu_confidence(asr_insitu)
  expect_equal(ncol(conf), 4)
})
