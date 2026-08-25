# Prepares a phylogeny for use in the insitu pipeline

Handles non-ultrametricity and incomplete taxon sampling. The tree must
be fully bifurcating.

## Usage

``` r
prep_phylo(phy, force_ultrametric = TRUE, sampling_fraction = 1)
```

## Arguments

- phy:

  A `phylo` object.

- force_ultrametric:

  Logical. If `TRUE` and the tree is not ultrametric, it is forced
  ultrametric using
  [`phytools::force.ultrametric`](https://rdrr.io/pkg/phytools/man/force.ultrametric.html).
  Default: `TRUE`.

- sampling_fraction:

  The proportion of extant species included in the phylogeny (between 0
  and 1). Used to correct branch length-based rate estimates. Default:
  `1` (complete sampling assumed).

## Value

The prepared `phylo` object, with `sampling_fraction` stored as an
attribute.
