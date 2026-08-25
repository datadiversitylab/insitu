# Flags and filters in-situ events relative to island formation age

Compares the age of each classified node to the formation age of the
island it is assigned to. Nodes that predate island formation are
biologically impossible and represent either errors in the ancestral
state reconstruction or lineages whose common ancestor actually predates
the island. These nodes are flagged and optionally removed. This
function also computes island-age-corrected rates by replacing the full
path length denominator with the time elapsed since island formation,
which is the relevant timescale for assessing how rapidly an island has
diversified.

## Usage

``` r
age_correct_events(events, phy, island_ages, remove_predating = FALSE)
```

## Arguments

- events:

  The data frame returned by `map_insitu_events`.

- phy:

  The matched phylogenetic tree.

- island_ages:

  A named numeric vector of island formation ages, in the same units as
  the branch lengths of `phy`. Names must match the island names in
  `events`.

- remove_predating:

  Logical. If `TRUE`, nodes that predate island formation are removed
  from the returned events table. If `FALSE` (default), they are
  retained but flagged in the `predates_island` column.

## Value

The events data frame with three additional columns:

- `node_age`:

  Age of the node in the same units as branch lengths (time before
  present).

- `island_age`:

  Formation age of the island assigned to that node.

- `predates_island`:

  Logical. `TRUE` if the node age exceeds the island formation age,
  meaning the event predates the island's existence.
