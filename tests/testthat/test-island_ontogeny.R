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

# Create named numeric vector of island ages
ages <- c(42,
          50,
          25,
          46,
          10,
          15)
island_names <- c("Cuba",
                  "Hispaniola",
                  "Jamaica",
                  "Puerto Rico",
                  "Cat Island",
                  "Ile de la Tortue")
# Assign the names to the numeric vector
names(ages) <- island_names

# For testing insitu_rate_by_island_age
new_events <- age_correct_events(events = map_events,
                                 phy = matched$phy,
                                 island_ages = ages,
                                 remove_predating = TRUE)

########

### age_correct_events ###

# Correct output class (data.frame)
test_that("age_correct_events returns a data.frame", {
  correct <- age_correct_events(events = map_events,
                                phy = matched$phy,
                                island_ages = ages,
                                remove_predating = FALSE)
  expect_s3_class(correct, "data.frame")
})

# Correct number of columns in the data.frame (9)
test_that("The resulting data.frame has the correct number of columns", {
  correct <- age_correct_events(events = map_events,
                                phy = matched$phy,
                                island_ages = ages,
                                remove_predating = FALSE)
  expect_equal(ncol(correct), 9)
})

# For this example data, two nodes should be flagged (remove_predating = FALSE)
test_that("The two nodes that meet the age criteria are flagged", {
  correct <- age_correct_events(events = map_events,
                                phy = matched$phy,
                                island_ages = ages,
                                remove_predating = FALSE)
  true_flags <- which(correct$predates_island == TRUE)
  expect_equal(length(true_flags), 2)
})

# For this example data, if remove_predating = TRUE, the number of rows should
#  be 22
test_that("When remove_predating = TRUE, the two flagged nodes should be removed", {
  correct <- age_correct_events(events = map_events,
                                phy = matched$phy,
                                island_ages = ages,
                                remove_predating = TRUE)
  expect_equal(nrow(correct), 22)
})

### insitu_rate_by_island_age ###

test_that("insitu_rate_by_island_age returns a data.frame", {
  rate_age <- insitu_rate_by_island_age(phy = matched$phy,
                                        events = new_events,
                                        island_ages = ages)
  expect_s3_class(rate_age, "data.frame")
})

test_that("The object returned by insitu_rate_by_island_age has a row per island", {
  rate_age <- insitu_rate_by_island_age(phy = matched$phy,
                                        events = new_events,
                                        island_ages = ages)
  expect_equal(nrow(rate_age), length(unique(matched$PAM$locale)))
})
