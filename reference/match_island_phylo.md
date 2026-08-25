# Imports and validates dated phylogenies with associated island occurrence data

Imports and validates dated phylogenies with associated island
occurrence data

## Usage

``` r
match_island_phylo(phy, locs, exclude = TRUE)
```

## Arguments

- phy:

  The phylogenetic tree associated with your data

- locs:

  A dataframe with 2 columns: "species" and "locale". Each row should
  represent a single locality where that species is located (e.x., if a
  species is found in multiple localities, there should be multiple rows
  for that species). This package focuses on island occurrences, so if a
  species is found on the mainland, write its locality as "Mainland".

- exclude:

  (Logical) Set to TRUE if you want to exclude mainland species from the
  output. Default: TRUE

## Value

A list with 3 items:

- `PAM`: A presence-absence matrix describing where species occur

- `sp_df`: A dataframe summarizing what species, and how many, are on
  each island

- `phy`: A pruned version of the user-provided phylogeny that includes
  only species that are in the user-provided locality data

## Examples

``` r
# Read tree trimmed from Patton et al. 2021
tree <- ape::read.tree(system.file("extdata",
                                   "Patton_etal_trimmed.tree",
                                    package = "insitu"))
# Read species location dataframe
dat <- read.csv(system.file("extdata",
                            "anolis_dat.csv",
                             package = "insitu"))
# Match island data and phylogeny
matched <- match_island_phylo(phy = tree, locs = dat)
#> ℹ Species dropped from the tree because they were not in the data: transversalis and heterodermus
#> ℹ Species dropped from the data because they were not in the tree: loysiana
#> ℹ Matching complete! Here are the first 5 locales in your data:
#> # A tibble: 5 × 3
#>   locale           species_list                                         richness
#>   <chr>            <chr>                                                   <int>
#> 1 Cat Island       distichus                                                   1
#> 2 Cuba             luteogularis, equestris, alutaceus, vanidicus, porc…       10
#> 3 Hispaniola       coelestinus, aliniger, bahorucoensis, olssoni, inso…       10
#> 4 Ile de la Tortue distichus                                                   1
#> 5 Jamaica          valencienni, lineatopus, garmani, grahami, sagrei           5
```
