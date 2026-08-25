# Warnings are triggered here due to the simplicity of the simulation parameters
#  in these test cases, so they have been suppressed.

test_that("The returned result is a data.frame", {
  test <- suppressWarnings(insitu_power(n_sim = 10, n_tips = 30,
                       birth_rate = 0.5, death_rate = 0,
                       colonization_rate = 0.3, island_extinction_fraction = 0,
                       n_islands = 3, threshold = 0.5, model = "ER",
                       use_simmap = FALSE, nsim = 2))

  expect_s3_class(test, "data.frame")
})

test_that("The returned result has 7 columns", {
  test <- suppressWarnings(insitu_power(n_sim = 10, n_tips = 30,
                                        birth_rate = 0.5, death_rate = 0,
                                        colonization_rate = 0.3, island_extinction_fraction = 0,
                                        n_islands = 3, threshold = 0.5, model = "ER",
                                        use_simmap = FALSE, nsim = 2))

  expect_equal(ncol(test), 7)
})
