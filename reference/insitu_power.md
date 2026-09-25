# Power analysis for in-situ speciation event recovery

Runs repeated simulations under known colonization, speciation, and
extinction parameters and measures how reliably the insitu pipeline
recovers the true in-situ events. For each simulation, a tree and PAM
are generated using
[`simulate_island`](https://datadiversitylab.github.io/insitu/reference/simulate_island.md),
the full pipeline is applied, and recovered events are compared to the
known ground truth.

## Usage

``` r
insitu_power(
  n_sim = 100,
  n_tips,
  diversification_rate,
  extinction_fraction = 0,
  diversification_rate_island = NULL,
  extinction_fraction_island = NULL,
  mainland_colonization_rate,
  inter_island_rate,
  n_islands = 1,
  n_mainland = 1,
  threshold = 0.5,
  model = "ER",
  use_simmap = FALSE,
  nsim = 10
)
```

## Arguments

- n_sim:

  Number of simulation replicates. Default: `100`.

- n_tips:

  Total number of tips per simulated tree.

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

- mainland_colonization_rate:

  Rate of mainland-to-island colonization events per unit of total
  mainland branch length.

- inter_island_rate:

  Rate of inter-island dispersal events per unit of total island branch
  length.

- n_islands:

  Number of islands to simulate. Default: `1`.

- n_mainland:

  Number of mainland tips. Default: `1`.

- threshold:

  ASR probability threshold passed to `map_insitu_events`. Default:
  `0.5`.

- model:

  Transition model passed to `run_geo_asr` or `simmap_insitu`. Default:
  `"ER"`.

- use_simmap:

  Logical. If `TRUE`, uses `simmap_insitu` instead of `run_geo_asr`.
  Default: `FALSE`.

- nsim:

  Number of stochastic maps per replicate when `use_simmap = TRUE`.
  Default: `10`.

## Value

A data frame with one row per simulation replicate and columns `sim`,
`n_true_insitu`, `n_recovered`, `n_false_neg`, `n_false_pos`,
`sensitivity`, and `precision`.
