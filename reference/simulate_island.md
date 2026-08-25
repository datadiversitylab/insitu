# Simulates island assemblages under known colonization, speciation, and extinction parameters

Generates a dated phylogeny and presence-absence matrix under a
birth-death-colonization process with known parameters, returning both
the simulated data and the true event table for use in
[`insitu_power`](https://datadiversitylab.github.io/insitu/reference/insitu_power.md).

## Usage

``` r
simulate_island(
  n_tips,
  birth_rate,
  death_rate = 0,
  colonization_rate,
  island_extinction_fraction = 0,
  n_islands = 1,
  seed = NULL
)
```

## Arguments

- n_tips:

  Total number of tips in the simulated tree.

- birth_rate:

  Per-lineage speciation rate.

- death_rate:

  Per-lineage extinction rate. Default: `0`.

- colonization_rate:

  Expected number of colonization events per unit of total branch
  length. Controls how many internal nodes are selected as colonization
  origins.

- island_extinction_fraction:

  Proportion of island tips to remove at random, simulating extinction.
  Default: `0`.

- n_islands:

  Number of distinct islands. Colonization nodes are distributed across
  islands in sequence. Default: `1`.

- seed:

  Optional random seed. Default: `NULL`.

## Value

A named list with elements `phy`, `PAM`, and `true_events`.
`true_events` has the same structure as the output of
`map_insitu_events` and serves as the ground truth for power analysis.

## Details

The simulation places colonization events directly on the tree by
randomly selecting internal nodes as colonization origins. The number of
colonization events is drawn from a Poisson distribution with mean equal
to `colonization_rate` multiplied by the total branch length of the
tree. Tips descending from a colonization node are assigned to an
island; all other tips are mainland. True in-situ events are nodes where
both descendant clades share exactly one island, computed directly from
the simulated topology without any ASR.
