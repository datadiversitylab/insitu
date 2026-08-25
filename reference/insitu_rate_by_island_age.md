# Computes island-age-corrected in-situ speciation rates

Modifies
[`insitu_speciation_rate`](https://datadiversitylab.github.io/insitu/reference/insitu_speciation_rate.md)
to use island formation age as the time denominator rather than total
root-to-tip path length. This gives the rate of in-situ speciation per
unit time that the island has existed. This approach allows for
comparing diversification dynamics across islands of different
geological ages.Events that predate island formation are excluded from
the count before computing rates.

## Usage

``` r
insitu_rate_by_island_age(phy, events, island_ages)
```

## Arguments

- phy:

  The matched phylogenetic tree.

- events:

  The data frame returned by
  [`age_correct_events`](https://datadiversitylab.github.io/insitu/reference/age_correct_events.md),
  which must include `node_age`, `island_age`, and `predates_island`
  columns.

- island_ages:

  A named numeric vector of island formation ages, in the same units as
  branch lengths.

## Value

A data frame with one row per island and columns:

- `island`:

  Island name.

- `n_insitu`:

  Number of valid in-situ events (excluding those predating the island).

- `island_age`:

  Island formation age.

- `insitu_rate_corrected`:

  In-situ events per unit time since island formation.
