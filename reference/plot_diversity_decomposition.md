# Stacked bar chart of in-situ, import, and export diversity per island

For each island, plots a stacked bar showing the number of in-situ
speciation, import, and export events as returned by
`decompose_diversity`.

## Usage

``` r
plot_diversity_decomposition(div)
```

## Arguments

- div:

  A data frame returned by `decompose_diversity`.

## Value

Invisibly returns the matrix used for plotting.
