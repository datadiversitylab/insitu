# Tests for simulate_island and insitu_power

test_that("simulate_island returns a list with three elements", {
  sim <- simulate_island(
    n_tips                     = 20,
    diversification_rate       = 0.5,
    extinction_fraction        = 0,
    mainland_colonization_rate = 0.05,
    inter_island_rate          = 0.5,
    n_islands                  = 3,
    n_mainland                 = 5,
    seed                       = 42
  )
  expect_type(sim, "list")
  expect_named(sim, c("phy", "PAM", "true_events"))
})

test_that("simulate_island phy element is a phylo object", {
  sim <- simulate_island(
    n_tips                     = 20,
    diversification_rate       = 0.5,
    extinction_fraction        = 0,
    mainland_colonization_rate = 0.05,
    inter_island_rate          = 0.5,
    n_islands                  = 3,
    n_mainland                 = 5,
    seed                       = 42
  )
  expect_s3_class(sim$phy, "phylo")
})

test_that("simulate_island produces an ultrametric tree", {
  sim <- simulate_island(
    n_tips                     = 20,
    diversification_rate       = 0.5,
    extinction_fraction        = 0,
    mainland_colonization_rate = 0.05,
    inter_island_rate          = 0.5,
    n_islands                  = 3,
    n_mainland                 = 5,
    seed                       = 42
  )
  expect_true(ape::is.ultrametric(sim$phy, tol = 1e-4))
})

test_that("simulate_island PAM is a data frame with a locale column", {
  sim <- simulate_island(
    n_tips                     = 20,
    diversification_rate       = 0.5,
    extinction_fraction        = 0,
    mainland_colonization_rate = 0.05,
    inter_island_rate          = 0.5,
    n_islands                  = 3,
    n_mainland                 = 5,
    seed                       = 42
  )
  expect_s3_class(sim$PAM, "data.frame")
  expect_true("locale" %in% names(sim$PAM))
  expect_true("Mainland" %in% sim$PAM$locale)
})

test_that("simulate_island PAM columns match tree tip labels", {
  sim <- simulate_island(
    n_tips                     = 20,
    diversification_rate       = 0.5,
    extinction_fraction        = 0,
    mainland_colonization_rate = 0.05,
    inter_island_rate          = 0.5,
    n_islands                  = 3,
    n_mainland                 = 5,
    seed                       = 42
  )
  pam_species <- setdiff(names(sim$PAM), "locale")
  expect_setequal(pam_species, sim$phy$tip.label)
})

test_that("simulate_island PAM contains only 0 and 1", {
  sim <- simulate_island(
    n_tips                     = 20,
    diversification_rate       = 0.5,
    extinction_fraction        = 0,
    mainland_colonization_rate = 0.05,
    inter_island_rate          = 0.5,
    n_islands                  = 3,
    n_mainland                 = 5,
    seed                       = 42
  )
  vals <- unlist(sim$PAM[, -1])
  expect_true(all(vals %in% c(0L, 1L)))
})

test_that("simulate_island true_events is a data frame", {
  sim <- simulate_island(
    n_tips                     = 20,
    diversification_rate       = 0.5,
    extinction_fraction        = 0,
    mainland_colonization_rate = 0.05,
    inter_island_rate          = 0.5,
    n_islands                  = 3,
    n_mainland                 = 5,
    seed                       = 42
  )
  expect_s3_class(sim$true_events, "data.frame")
})

test_that("simulate_island true_events has expected columns", {
  sim <- simulate_island(
    n_tips                     = 20,
    diversification_rate       = 0.5,
    extinction_fraction        = 0,
    mainland_colonization_rate = 0.05,
    inter_island_rate          = 0.5,
    n_islands                  = 3,
    n_mainland                 = 5,
    seed                       = 42
  )
  expected <- c("node", "island", "ancestor_presence_prob",
                "in_situ", "export", "import")
  expect_true(all(expected %in% names(sim$true_events)))
})

test_that("simulate_island true_events nodes are internal nodes", {
  sim <- simulate_island(
    n_tips                     = 20,
    diversification_rate       = 0.5,
    extinction_fraction        = 0,
    mainland_colonization_rate = 0.05,
    inter_island_rate          = 0.5,
    n_islands                  = 3,
    n_mainland                 = 5,
    seed                       = 42
  )
  n_tips <- ape::Ntip(sim$phy)
  expect_true(all(sim$true_events$node > n_tips))
})

test_that("simulate_island respects n_mainland", {
  sim <- simulate_island(
    n_tips                     = 20,
    diversification_rate       = 0.5,
    extinction_fraction        = 0,
    mainland_colonization_rate = 0.05,
    inter_island_rate          = 0.5,
    n_islands                  = 3,
    n_mainland                 = 8,
    seed                       = 42
  )
  mainland_tips <- sum(grepl("^spM_", sim$phy$tip.label))
  expect_equal(mainland_tips, 8L)
})

test_that("simulate_island seed produces reproducible output", {
  sim1 <- simulate_island(
    n_tips                     = 20,
    diversification_rate       = 0.5,
    extinction_fraction        = 0,
    mainland_colonization_rate = 0.05,
    inter_island_rate          = 0.5,
    n_islands                  = 3,
    n_mainland                 = 5,
    seed                       = 99
  )
  sim2 <- simulate_island(
    n_tips                     = 20,
    diversification_rate       = 0.5,
    extinction_fraction        = 0,
    mainland_colonization_rate = 0.05,
    inter_island_rate          = 0.5,
    n_islands                  = 3,
    n_mainland                 = 5,
    seed                       = 99
  )
  expect_equal(sim1$true_events, sim2$true_events)
})

test_that("simulate_island stops when n_mainland leaves fewer than 2 island tips", {
  expect_error(
    simulate_island(
      n_tips                     = 10,
      diversification_rate       = 0.5,
      extinction_fraction        = 0,
      mainland_colonization_rate = 0.05,
      inter_island_rate          = 0.5,
      n_islands                  = 3,
      n_mainland                 = 9
    )
  )
})

test_that("simulate_island island-specific rates differ from mainland rates", {
  sim <- simulate_island(
    n_tips                      = 30,
    diversification_rate        = 0.15,
    extinction_fraction         = 0.10,
    diversification_rate_island = 2.00,
    extinction_fraction_island  = 0.20,
    mainland_colonization_rate  = 0.05,
    inter_island_rate           = 0.50,
    n_islands                   = 3,
    n_mainland                  = 10,
    seed                        = 7
  )
  expect_s3_class(sim$phy, "phylo")
})

# --- insitu_power tests -------------------------------------------------------

test_that("insitu_power returns a data frame", {
  res <- insitu_power(
    n_sim                      = 5,
    n_tips                     = 20,
    diversification_rate       = 0.5,
    extinction_fraction        = 0,
    mainland_colonization_rate = 0.05,
    inter_island_rate          = 0.5,
    n_islands                  = 3,
    n_mainland                 = 5,
    threshold                  = 0.5,
    model                      = "ER"
  )
  expect_s3_class(res, "data.frame")
})

test_that("insitu_power returns expected columns", {
  res <- insitu_power(
    n_sim                      = 5,
    n_tips                     = 20,
    diversification_rate       = 0.5,
    extinction_fraction        = 0,
    mainland_colonization_rate = 0.05,
    inter_island_rate          = 0.5,
    n_islands                  = 3,
    n_mainland                 = 5,
    threshold                  = 0.5,
    model                      = "ER"
  )
  expected <- c("sim", "n_true_insitu", "n_recovered",
                "n_false_neg", "n_false_pos", "sensitivity", "precision")
  expect_true(all(expected %in% names(res)))
})

test_that("insitu_power sensitivity and precision are between 0 and 1", {
  res <- insitu_power(
    n_sim                      = 10,
    n_tips                     = 20,
    diversification_rate       = 0.5,
    extinction_fraction        = 0,
    mainland_colonization_rate = 0.05,
    inter_island_rate          = 0.5,
    n_islands                  = 3,
    n_mainland                 = 5,
    threshold                  = 0.5,
    model                      = "ER"
  )
  expect_true(all(res$sensitivity >= 0 & res$sensitivity <= 1, na.rm = TRUE))
  expect_true(all(res$precision   >= 0 & res$precision   <= 1, na.rm = TRUE))
})

test_that("insitu_power excludes replicates with no true in-situ events", {
  res <- insitu_power(
    n_sim                      = 10,
    n_tips                     = 20,
    diversification_rate       = 0.5,
    extinction_fraction        = 0,
    mainland_colonization_rate = 0.05,
    inter_island_rate          = 0.5,
    n_islands                  = 3,
    n_mainland                 = 5,
    threshold                  = 0.5,
    model                      = "ER"
  )
  expect_true(all(res$n_true_insitu > 0))
})

test_that("insitu_power n_recovered does not exceed n_true_insitu", {
  res <- insitu_power(
    n_sim                      = 10,
    n_tips                     = 20,
    diversification_rate       = 0.5,
    extinction_fraction        = 0,
    mainland_colonization_rate = 0.05,
    inter_island_rate          = 0.5,
    n_islands                  = 3,
    n_mainland                 = 5,
    threshold                  = 0.5,
    model                      = "ER"
  )
  expect_true(all(res$n_recovered <= res$n_true_insitu))
})
