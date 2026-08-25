# Statistical comparison of diversity decomposition metrics across multiple island systems

Tests whether in-situ speciation, import, and export counts differ
across islands using a chi-square test, and whether in-situ and
colonization rates differ across islands using a Kruskal-Wallis test.
Rates are computed internally from the counts and total species per
island.

## Usage

``` r
compare_islands(div)
```

## Arguments

- div:

  A data frame returned by `decompose_diversity`.

## Value

A named list with elements:

- `counts_test`:

  Chi-square test comparing n_insitu, n_import, and n_export across
  islands.

- `insitu_rate_test`:

  Kruskal-Wallis test on per-island in-situ speciation rates (n_insitu /
  n_species).

- `colonization_rate_test`:

  Kruskal-Wallis test on per-island colonization rates ((n_import +
  n_export) / 2 / n_species).
