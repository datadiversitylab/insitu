#' Imports and validates dated phylogenies with associated island occurrence data
#'
#'
#' @param phy The phylogenetic tree associated with your data
#' @param locs A dataframe with 2 columns: "species" and "locale". Each row
#' should represent a single locality where that species is located (e.x., if
#' a species is found in multiple localities, there should be multiple rows for
#' that species). This package focuses on island occurrences, so if a species is
#' found on the mainland, write its locality as "Mainland".
#' @param exclude (Logical) Set to TRUE if you want to exclude mainland species
#' from the output. Default: TRUE
#' @return A list with 3 items:
#' * `PAM`: A presence-absence matrix describing where species occur
#' * `sp_df`: A dataframe summarizing what species, and how many, are on each island
#' * `phy`: A pruned version of the user-provided phylogeny that includes only species that are in the user-provided locality data
#' @examples
#' # Read tree trimmed from Patton et al. 2021
#' tree <- ape::read.tree(system.file("extdata",
#'                                    "Patton_etal_trimmed.tree",
#'                                     package = "insitu"))
#' # Read species location dataframe
#' dat <- read.csv(system.file("extdata",
#'                             "anolis_dat.csv",
#'                              package = "insitu"))
#' # Match island data and phylogeny
#' matched <- match_island_phylo(phy = tree, locs = dat)
#'
#' @export

match_island_phylo <- function(phy, locs, exclude = TRUE){
  # If the user wants to exclude mainland species, do so
  if(exclude == TRUE){
    locs <- locs[which(locs$locale != "Mainland"),]
  }

  # Make sure that each species in locs is also present in the tree
  # If not, remove the species that don't match (and make a note about it)

  # Create a dataframe with just species names for ease of using treedata
  temp_df <- as.data.frame(unique(locs$species))
  rownames(temp_df) <- temp_df[,1]

  td <- geiger::treedata(phy, temp_df, warnings = FALSE)

  # Write a nice message about which species were dropped in each case
  # Difference between data and phy
  # (the tips that were dropped from phy because they were not in data)
  drop_sp <- phy$tip.label[!(phy$tip.label %in% td$phy$tip.label)]
  # Write the message
  if (!getOption("insitu.silent", FALSE)) {
    cli::cli_alert_info("Species dropped from the tree because they were not in the data: {drop_sp}")
  }

  # Difference between phy and data
  # (species that were dropped from data because they were not in phy)
  # Dataframe from treedata (with dropped species)
  temp_dat <- as.data.frame(td$data)
  colnames(temp_dat) <- "species"

  drop_sp <- temp_df[,1][!(temp_df[,1] %in% temp_dat$species)]
  # Write the message
  if (!getOption("insitu.silent", FALSE)) {
    cli::cli_alert_info("Species dropped from the data because they were not in the tree: {drop_sp}")
  }

  # Save the pruned tree
  phy <- td$phy

  # Remove species that were not in the tree from the locality data
  final_locs <- dplyr::inner_join(x = locs, y = temp_dat, by = "species")

  # Next, create the PAM and sp_df objects to return
  # Create PAM
  PAM <- reshape2::dcast(
    final_locs,
    locale ~ species,
    fun.aggregate = length,
    value.var = "species"
  )

  # Create sp_df, with a list of all species on each island and a count
  sp_df <- final_locs |>
    # First, group by island
    dplyr::group_by(.data$locale) |>
    # Next, paste all of the species into one list and count them
    dplyr::summarize(
      species_list = paste(.data$species, collapse = ", "),
      richness = dplyr::n_distinct(.data$species)
    )

  results <- list(
    PAM = PAM,
    sp_df = sp_df,
    phy = phy
  )

  # Print out head of sp_df to allow the user to quickly check it
  if (!getOption("insitu.silent", FALSE)) {
    cli::cli_alert_info("Matching complete! Here are the first 5 locales in your data:")
    print(utils::head(sp_df, n = 5))
  }

  return(results)
}
