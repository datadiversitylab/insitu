# Simulates island assemblages under known colonization, speciation, and extinction parameters

Generates a dated phylogeny and presence-absence matrix under a
birth-death-colonization process with known parameters, returning both
the simulated data and the true event table for use in
[`insitu_power`](https://datadiversitylab.github.io/insitu/reference/insitu_power.md).

## Usage

``` r
simulate_island(
  n_tips,
  diversification_rate,
  extinction_fraction = 0,
  diversification_rate_island = NULL,
  extinction_fraction_island = NULL,
  colonization_rate,
  n_islands = 1,
  n_mainland = 1,
  monophyletic_island = TRUE,
  seed = NULL
)
```

## Arguments

- n_tips:

  Total number of tips in the simulated tree.

- diversification_rate:

  Net diversification rate (birth - death).

- extinction_fraction:

  Relative extinction rate (death / birth), between 0 and 1. Default:
  `0`.

- diversification_rate_island:

  Island-specific net diversification rate. Defaults to
  `diversification_rate` when not set.

- extinction_fraction_island:

  Island-specific relative extinction. Defaults to `extinction_fraction`
  when not set.

- colonization_rate:

  Expected number of inter-island colonization events per unit of total
  branch length.

- n_islands:

  Number of distinct islands. Default: `1`.

- n_mainland:

  Number of tips forced to be mainland. Default: `1`.

- monophyletic_island:

  Logical. If `TRUE` (default), all island species descend from a single
  crown node. If `FALSE`, colonization nodes are placed freely across
  the tree.

- seed:

  Optional random seed. Default: `NULL`.

## Value

A named list with elements `phy`, `PAM`, and `true_events`.

## Details

When `monophyletic_island = TRUE`, all island species descend from a
single crown node, reflecting the biology of most island radiations
where all island species trace back to a single colonization from the
mainland. Inter-island colonization events are then placed within that
clade. When `monophyletic_island = FALSE`, colonization nodes are placed
freely across the tree.

Speciation and extinction rates are specified as net diversification
rate and extinction fraction. Birth and death rates are derived
internally as: `birth = diversification / (1 - extinction_fraction)` and
`death = birth * extinction_fraction`.
