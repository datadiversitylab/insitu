#' Assigns tips as island endemic or non-endemic based on locality data
#'
#' This function classifies species both on the island-level and on
#' the archipelago-level. For example, a species that occurs in two islands within
#' the same archipelago is not endemic on an island-level (0 in the island column),
#' but is endemic on an archipelago-level (1 in the archipelago column).
#'
#' @param occs Occurrence data (in a data.frame) that includes information about
#' several species and where they occur. This function assumes that there is a column
#' called "species" where the species name is written, and a column called "locale"
#' where the corresponding locale is written. If one species occupies multiple locales,
#' there should be multiple rows for that species in this dataframe.
#' @param islands A vector of island names that help classify each tip.
#'
#'
#' @export

classify_tips <- function(occs, islands){
  # Check the occs dataframe and determine:
  # 1. How many locales each species occupies
  # 2. Whether each species is island endemic or non-endemic at each level
  #    (island vs archipelago)

  # How many locales does each species occupy?
  occs_locality <- occs |>
    # Group by species to account for one species having multiple rows
    dplyr::group_by(.data$species) |>
    # Add all localities for each species, with a comma separating them
    dplyr::summarize(locality = paste(.data$locale, collapse = ", "))

  # Create a new dataframe that will be populated with info about endemicity
  endemic <- as.data.frame(matrix(nrow = length(occs_locality$species),
                                  ncol = 3))
  # Add the species column as the first column
  endemic[,1] <- occs_locality$species
  # Change column names
  colnames(endemic) <- c("species", "archipelago", "island")

  # For each species, determine whether it's endemic at each level
  #  (island and archipelago)
  for(i in c(1:length(occs_locality$locality))){
    # Split the string by commas, and remove extra spacing
    locale_vec <- stringr::str_split(occs_locality$locality[i], ",\\s*")[[1]]

    # Check if all of the islands are within the archipelago
    if(all(locale_vec %in% islands)){
      endemic$archipelago[i] <- 1
    } else {
      endemic$archipelago[i] <- 0
    }

    # If locale_vec is more than one entry, then the species is not endemic at
    #  the island level
    if(length(locale_vec) > 0){
      endemic$island[i] <- 0
    } else {
      endemic$island[i] <- 1
    }
  } # End for
}
