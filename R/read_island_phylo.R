#' Imports and validates dated phylogenies with associated island occurrence data
#'
#' Requirements: make sure that the row names for your data are the species names
#'
#' @param phy_path The local file path to the tree file
#' @param occ_path The local file path to the occurrence data file
#' @return A list with 2 items: the pruned phylogenetic tree and associated data
#'
#' @export

read_island_phylo <- function(phy_path, occ_path){
  # Read in tree
  phy <- ape::read.tree(phy_path)

  # Read in occurrence data
  dat <- utils::read.csv(occ_path, row.names = 1)

  # Check that the tips and data match
  suppressWarnings(td <- geiger::treedata(phy = phy, data = dat))

  # Return objects, with dat as a data.frame
  new_dat <- as.data.frame(td$data)

  return(list(phy = td$phy,
              dat = new_dat))
}
