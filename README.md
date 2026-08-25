# insitu
An R package for detecting *in situ* speciation and providing a clearer picture of how biodiversity emerged on isolated systems.

## Overview
A speciation event is considered to be *in situ* when it occurs within a specific geographical range, as opposed to a lineage colonizing a location from elsewhere. The importance of *in situ* speciation for understanding diversity patterns has been emphasized by decades of work, much of which discusses how *in situ* speciation contributes to species richness on islands.

The `insitu` R package provides a suite of functions that allows the user to estimate *in situ* speciation events and quantify their contribution to biodiversity in discrete geographical units, with a primary focus on island systems and island-like environments. This vignette will cover the functions from the `insitu` package that are necessary to perform a basic analysis that classifies nodes on a phylogeny as *in situ* speciation, colonization, or dispersal events, as well as computing diversity decomposition metrics and visualizing results.

## Installation
The `insitu` R package can be installed using `remotes`:
```r
remotes::install_github("datadiversitylab/insitu")
```

## Getting Started
Please find our vignette describing the main functions of `insitu` and how to get started using it for your analyses [on our `pkgdown` website linked here](https://datadiversitylab.github.io/insitu/articles/GettingStarted.html).
