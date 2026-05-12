# Read in example data and tree
# Read tree trimmed from Patton et al. 2021
tree <- ape::read.tree(system.file("extdata",
                                   "Patton_etal_trimmed.tree",
                                    package = "insitu"))

# Read species location dataframe
dat <- read.csv(system.file("extdata",
                            "anolis_dat.csv",
                             package = "insitu"))

########

# Checking types
# Full list
test_that("match_island_phylo returns a list", {
  matched <- match_island_phylo(phy = tree, locs = dat)
  expect_type(matched, "list")
})

# First element (data.frame)
test_that("The first element of the returned list is a data.frame", {
  matched <- match_island_phylo(phy = tree, locs = dat)
  expect_s3_class(matched$PAM, "data.frame")
})

# Second element (data.frame)
test_that("The second element of the returned list is a data.frame", {
  matched <- match_island_phylo(phy = tree, locs = dat)
  expect_s3_class(matched$sp_df, "data.frame")
})

# Third element (list)
test_that("The third element of the returned list is of type list", {
  matched <- match_island_phylo(phy = tree, locs = dat)
  expect_type(matched$phy, "list")
})
