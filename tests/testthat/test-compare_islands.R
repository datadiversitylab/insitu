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

# Decompose diversity
div <- decompose_diversity(phy = matched$phy,
                           events = map_events,
                           PAM = matched$PAM)

test <- compare_islands(div)

########
# Function returns a list
test_that("compare_islands returns a list", {
  # This function throws a chi-squared warning that is not indiciative of
  #  anything actually breaking
  comp <- suppressWarnings(compare_islands(div))
  expect_type(comp, "list")
})

# The first element of the list has length 9
test_that("The first list element (counts_test) has length 9", {
  comp <- suppressWarnings(compare_islands(div))
  expect_equal(length(comp[[1]]), 9)
})

# The second element of the list has length 5
test_that("The first list element (insitu_rate_test) has length 5", {
  comp <- suppressWarnings(compare_islands(div))
  expect_equal(length(comp[[2]]), 5)
})

# The third element of the list has length 5
test_that("The first list element (colonization_rate_test) has length 5", {
  comp <- suppressWarnings(compare_islands(div))
  expect_equal(length(comp[[3]]), 5)
})
