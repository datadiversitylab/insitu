test_that("simulate_island returns a list", {
  test <- simulate_island(n_tips = 10,
                          birth_rate = 0.5,
                          death_rate = 0,
                          colonization_rate = 0.3,
                          island_extinction_fraction = 0,
                          n_islands = 3,
                          seed = 42)
  expect_type(test, "list")
})

test_that("The first element of the returned list is of type list", {
  test <- simulate_island(n_tips = 10,
                          birth_rate = 0.5,
                          death_rate = 0,
                          colonization_rate = 0.3,
                          island_extinction_fraction = 0,
                          n_islands = 3,
                          seed = 42)
  expect_type(test$phy, "list")
})

test_that("The second element of the returned list is a data.frame", {
  test <- simulate_island(n_tips = 10,
                          birth_rate = 0.5,
                          death_rate = 0,
                          colonization_rate = 0.3,
                          island_extinction_fraction = 0,
                          n_islands = 3,
                          seed = 42)
  expect_s3_class(test$PAM, "data.frame")
})

test_that("The third element of the returned list is a data.frame", {
  test <- simulate_island(n_tips = 10,
                          birth_rate = 0.5,
                          death_rate = 0,
                          colonization_rate = 0.3,
                          island_extinction_fraction = 0,
                          n_islands = 3,
                          seed = 42)
  expect_s3_class(test$true_events, "data.frame")
})

