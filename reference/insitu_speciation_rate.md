# Estimates the rate of in-situ speciation events per lineage per unit time

For each tip, counts the number of in-situ speciation events on the path
from the root to the tip and divides that by the total branch length of
the path. When the tree was prepared with
[`prep_phylo`](https://datadiversitylab.github.io/insitu/reference/prep_phylo.md)
and a `sampling_fraction` less than 1 was set, branch lengths are scaled
upward by dividing by the sampling fraction before computing rates. This
approach corrects for the fact that incomplete sampling shortens the
total branch length of the observed tree relative to the true underlying
tree.

## Usage

``` r
insitu_speciation_rate(phy, events)
```

## Arguments

- phy:

  The matched phylogenetic tree prepared with
  [`prep_phylo`](https://datadiversitylab.github.io/insitu/reference/prep_phylo.md)

- events:

  The data frame returned by `map_insitu_events`.

## Value

A data frame with one row per tip and columns:

- `species`:

  Tip label.

- `n_insitu`:

  Number of in-situ speciation events on the path from root to tip.

- `path_length`:

  Total branch length from root to tip.

- `insitu_rate`:

  In-situ speciation events per unit time.
