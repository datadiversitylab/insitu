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
  mainland_colonization_rate,
  inter_island_rate,
  n_islands = 1,
  n_mainland = 1,
  seed = NULL
)
```

## Arguments

- n_tips:

  Total number of tips in the simulated tree.

- diversification_rate:

  Net mainland diversification rate (birth - death).

- extinction_fraction:

  Mainland relative extinction (death / birth), between 0 and 1.
  Default: `0`.

- diversification_rate_island:

  Island net diversification rate. Defaults to `diversification_rate`
  when not set.

- extinction_fraction_island:

  Island relative extinction. Defaults to `extinction_fraction` when not
  set.

- mainland_colonization_rate:

  Rate of mainland-to-island colonization events per unit of total
  mainland branch length.

- inter_island_rate:

  Rate of inter-island dispersal events per unit of total island branch
  length.

- n_islands:

  Number of distinct islands. Default: `1`.

- n_mainland:

  Number of mainland tips. Default: `1`.

- seed:

  Optional random seed. Default: `NULL`.

## Value

A named list with elements `phy`, `PAM`, and `true_events`.

## Details

The simulation models three distinct processes. Mainland lineages
colonize the island system at rate `mainland_colonization_rate`: the
number of independent colonization events is drawn from a Poisson
distribution with mean equal to that rate multiplied by the total
mainland branch length. When only one event is drawn, all island species
are monophyletic; when multiple events are drawn, the island assemblage
is polyphyletic. Once on the island system, lineages disperse between
islands at rate `inter_island_rate`. Speciation and extinction proceed
at island-specific rates that may differ from the mainland.

True in-situ events are nodes where both descendant clades share exactly
one island AND descend from the same mainland colonization event. Nodes
spanning two independent colonization-derived clades are not in-situ at
the system level even if their descendants share an island.
