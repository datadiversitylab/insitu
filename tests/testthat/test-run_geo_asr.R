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

########

test_that("run_geo_asr returns a list", {
  res <- run_geo_asr(phy = matched$phy, PAM = matched$PAM, model = "ER")
  expect_type(res, "list")
})

# Here, the list has 6 elements
test_that("The resulting list has 6 elements", {
  res <- run_geo_asr(phy = matched$phy, PAM = matched$PAM, model = "ER")
  expect_equal(length(res), 6)
})
