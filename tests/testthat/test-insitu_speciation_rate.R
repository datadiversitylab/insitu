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
test_that("insitu_speciation_rate returns a data.frame", {
  rates <- insitu_speciation_rate(phy = matched$phy, events = map_events)
  expect_s3_class(rates, "data.frame")
})

# Correct number of rows in the data.frame (should equal the number of tips)
test_that("The resulting data.frame has the correct number of rows", {
  rates <- insitu_speciation_rate(phy = matched$phy, events = map_events)
  expect_equal(nrow(rates), length(matched$phy$tip.label))
})
