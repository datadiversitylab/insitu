set.seed(42)

# Generate a random, non-ultrametric tree
nu_phy <- ape::rtree(n = 30)

# Read in ultrametric Anolis tree
u_phy <- ape::read.tree(system.file("extdata/Patton_etal_trimmed.tree",
                                    package = "insitu"))
########

test_that("prep_phylo works with ultrametric tree and reasonable sampling fraction", {
  expect_no_error(prep_phylo(u_phy))
})

test_that("prep_phylo issues a message when the tree is non-ultrametric (force_ultrametric = TRUE)", {
  expect_message(prep_phylo(nu_phy))
})

test_that("prep_phylo issues a warning when the tree is non-ultrametric (force_ultrametric = FALSE)", {
  expect_warning(prep_phylo(nu_phy, force_ultrametric = FALSE), "Tree is not ultrametric. Some functions in insitu may fail. Set force_ultrametric = TRUE to resolve this automatically.")
})

test_that("prep_phylo issues an error when the sampling fraction is illogical", {
  expect_error(prep_phylo(u_phy, sampling_fraction = 2))
})

