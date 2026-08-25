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

# Normally, you would want to run the function in the test cases themselves, but
#  because simmap_insitu takes a while to run, let's just do that once
test <- simmap_insitu(phy = matched$phy, PAM = matched$PAM, model = "ER", nsim = 2)

########

test_that("The result from simmap_insitu is a list", {
  expect_type(test, "list")
})

test_that("The resulting list from simmap_insitu is of length 6", {
  expect_equal(length(test), 6)
})
