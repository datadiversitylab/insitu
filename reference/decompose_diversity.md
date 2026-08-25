# Partition island species richness into in-situ, import, and export components

For each island, counts in-situ speciation, import, and export events.
In-situ speciation events are identified as cases in which an ancestral
node is assigned to the same region as its descendant. Import events are
cases where a descendant node is assigned to a different region than its
ancestor. Export events are the reverse: the ancestral node is assigned
to the focal island while its descendant is assigned to a different
region. Both ancestor and descendant must be assigned to a region for
any event to be detected.

## Usage

``` r
decompose_diversity(phy, events, PAM)
```

## Arguments

- phy:

  The matched phylogenetic tree.

- events:

  The data frame returned by `map_insitu_events`, with columns `node`,
  `island`, `in_situ`, `export`, and `import`.

- PAM:

  The presence-absence matrix with islands as rows and species as
  columns, with a `locale` column for island names.

## Value

A data frame with one row per island and columns:

- `island`:

  Island name.

- `n_species`:

  Total extant species on the island.

- `n_insitu`:

  Species derived from in-situ speciation.

- `n_import`:

  Import events into this island.

- `n_export`:

  Export events out of this island.

- `prop_insitu`:

  Proportion of species from in-situ speciation.
