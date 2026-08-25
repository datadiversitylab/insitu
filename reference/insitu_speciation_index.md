# In-situ speciation index per island

Computes a standardized index of in-situ contribution to diversity for
each island. The index (isci) is the number of in-situ speciation events
on an island divided by the total branch length of the island's subtree.
This is, the sum of all branch lengths of all descendant tips that are
exclusively assigned to that island.

## Usage

``` r
insitu_speciation_index(phy, events)
```

## Arguments

- phy:

  Matched phylogenetic tree

- events:

  The data frame returned by `map_insitu_events`

## Value

A data frame with one row per island and columns:

- `island`:

  Island name.

- `n_insitu`:

  Number of in-situ speciation events.

- `island_bl`:

  Total branch length of the island subtree.

- `isci`:

  In-situ contribution index: `n_insitu / island_bl`.
