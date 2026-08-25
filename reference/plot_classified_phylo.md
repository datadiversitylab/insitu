# Dated phylogeny with branches and nodes colored by in-situ vs. colonization classification

Dated phylogeny with branches and nodes colored by in-situ vs.
colonization classification

## Usage

``` r
plot_classified_phylo(phy, events, PAM)
```

## Arguments

- phy:

  The phylogenetic tree associated with your data

- events:

  The dataframe returned by `map_insitu_events`. If the user would like
  to input a custom dataframe, ensure that it has the following columns:
  node (internal node numbers associated with the given phylogenetic
  tree), island (the name of the island associated with the ancestral
  state for that node), and in_situ (logical - whether in situ
  speciation occurred at the given node).

- PAM:

  The presence-absence matrix associated with your data as returned by
  the `match_island_phylo` function.
