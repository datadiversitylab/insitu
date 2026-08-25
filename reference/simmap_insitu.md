# Uses stochastic character mapping to propagate ASR uncertainty through speciation classifications. Posterior distribution of in-situ event counts

When multiple trees are provided, the function propagates ASR
uncertainty from the stochastic mapping itself, and phylogenetic
uncertainty from the distribution of trees.

## Usage

``` r
simmap_insitu(phy, PAM, model = "ER", nsim = 10)
```

## Arguments

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

- model:

  The transition model to use in
  [`phytools::make.simmap`](https://rdrr.io/pkg/phytools/man/make.simmap.html).
  See documentation for
  [`ape::ace`](https://rdrr.io/pkg/ape/man/ace.html) for more
  information. Default: ER

- nsim:

  The number of simulations to use in
  [`phytools::make.simmap`](https://rdrr.io/pkg/phytools/man/make.simmap.html)
  Default: 10

## Value

A list of island-specific ancestral state reconstructions
