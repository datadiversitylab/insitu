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
  birth_rate,
  death_rate = 0,
  colonization_rate,
  island_extinction_fraction = 0,
  n_islands = 1,
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

- birth_rate:

  Per-lineage speciation rate.

- death_rate:

  Per-lineage extinction rate. Default: `0`.

- colonization_rate:

  Rate of mainland-to-island transitions per unit branch length.

- island_extinction_fraction:

  Proportion of island tips to remove, simulating extinction. Default:
  `0`.

- n_islands:

  Number of islands to simulate. Default: `1`.

- threshold:

  ASR probability threshold passed to `map_insitu_events`. Default:
  `0.5`.

- model:

  Transition model passed to `run_geo_asr` or `simmap_insitu`. Default:
  `"ER"`.

- use_simmap:

  Logical. If `TRUE`, uses `simmap_insitu` instead of `run_geo_asr` to
  propagate ASR uncertainty through the power analysis. Default:
  `FALSE`.

- nsim:

  Number of stochastic maps per replicate when `use_simmap = TRUE`.
  Default: `10`.

## Value

A data frame with one row per simulation replicate and columns:

- `sim`:

  Replicate number.

- `n_true_insitu`:

  Number of true in-situ events in the simulated data.

- `n_recovered`:

  Number of true in-situ events correctly recovered by the pipeline.

- `n_false_positive`:

  Nodes classified as in-situ that were not true in-situ events.

- `sensitivity`:

  True positive rate: `n_recovered / n_true_insitu`.

- `false_positive_rate`:

  False positive rate: ` n_false_positive / n_non_insitu_nodes`.

## Details

Sensitivity measures the proportion of true in-situ events that the
pipeline correctly identifies. Specificity measures the proportion of
true non-in-situ nodes that are correctly left unclassified as in-situ.
Both vary with tree size, extinction rate, and the ASR threshold, so
running this function across a range of parameter values directly
answers how trustworthy a given set of in-situ classifications is likely
to be.
