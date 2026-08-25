# Applies decision rules to reconstructed nodes to classify each speciation event as in-situ, colonization, or ambiguous

Applies decision rules to reconstructed nodes to classify each
speciation event as in-situ, colonization, or ambiguous

## Usage

``` r
map_insitu_events(recons, phy, PAM, threshold = 0.5)
```

## Arguments

- recons:

  A list of ancestral state reconstructions returned by `run_geo_asr`.
  We do not recommend inputting a custom list for this function, but if
  you need to, make sure that the list is composed of
  [`ape::ace`](https://rdrr.io/pkg/ape/man/ace.html) objects, each named
  after a locality. Refer to `run_geo_asr` for more information about
  how the original function estimates ancestral states.

- phy:

  The phylogenetic tree associated with your data

- PAM:

  A presence-absence matrix reflecting where each species of interest is
  located. The `match_island_phylo` function will create this PAM, but
  if the user would prefer to input a custom PAM, make sure that it has
  a column titled "locale" with the name of an island in each row, and
  each subsequent column is titled with a species name. Each species
  column should include either a 0 (absence) or a 1 (presence),
  signifying whether that species occurs on the island in a given row.

- threshold:

  The threshold of probability at which it is reasonable to infer that a
  given ancestral node occurred on the same island as a given set of
  tips. Default: 0.5

## Value

A dataframe including node numbers, islands associated with those nodes,
and whether an in situ speciation event was likely at that node/island
combination.
