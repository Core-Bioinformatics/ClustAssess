# Computes the NN adjacency matrix given the neighbours

Computes the NN adjacency matrix given the neighbours

## Usage

``` r
getNNmatrix(nnRanked, k = -1L, start = 0L, prune = 0)
```

## Arguments

- nnRanked:

  A matrix with the lists of the nearest neighbours for each point

- k:

  The number of neighbours to consider. Defaults to `-1`, which means
  all neighbours.

- start:

  The index of the first neighbour to consider. Defaults to `0`.

- prune:

  The threshold to prune the SNN matrix. If -1, the function will only
  return the NN matrix. Defaults to `0`.

## Value

A list with the NN and SNN adjacency matrices.
