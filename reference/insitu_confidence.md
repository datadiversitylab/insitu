# Computes per-node probabilities of island presence for all internal nodes

Computes per-node probabilities of island presence for all internal
nodes

## Usage

``` r
insitu_confidence(recons)
```

## Arguments

- recons:

  A list of ancestral state reconstructions returned by `run_geo_asr`.

## Value

A data frame with one row per node per island and columns:

- `node`:

  Internal node number.

- `island`:

  Island name.

- `prob`:

  Probability of the ancestor being present on the island.

- `support`:

  Qualitative support level: `"high"` (prob \>= 0.75), `"moderate"`
  (prob \>= 0.5), or `"low"` (prob \< 0.5).
