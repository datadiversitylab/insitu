# Estimates colonization rates events per lineage per unit time

For each tip, counts the number of dispersal events on the path from the
root to the tip and divides by the total branch length of that path.
Each dispersal event is recorded as both an import and an export, so the
total transitions divided by two to avoid double-counting. Accounts for
incomplete sampling fraction when available.

## Usage

``` r
colonization_rate(phy, events)
```

## Arguments

- phy:

  The marched phylogenetic tree.

- events:

  The data frame returned by `map_insitu_events`.

## Value

A data frame with one row per tip and columns:

- `species`:

  Tip label.

- `n_import`:

  Number of import events on the path from root to tip.

- `n_export`:

  Number of export events on the path from root to tip.

- `n_colonization`:

  Number of dispersal events

- `path_length`:

  Total branch length from root to tip.

- `colonization_rate`:

  Colonization events per unit time.
