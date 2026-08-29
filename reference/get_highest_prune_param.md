# Calculate the highest pruning parameter for the SNN graph given NN matrix

Given a NN adjacency matrix, the function calculates the highest pruning
parameter for the SNN graph that preserves the connectivity of the
graph.

## Usage

``` r
get_highest_prune_param(nn_matrix, n_neigh)
```

## Arguments

- nn_matrix:

  The adjacency matrix of the nearest neighbour graph.

- n_neigh:

  The number of nearest neighbours.

## Value

A list with the following fields:

- `prune_value`: The value of the highest pruning parameter.

- `adj_matrix`: The adjacency matrix of the SNN graph after pruning.

## Note

Given the way the SNN graph is built, the possible values for the
pruning parameter are limited and can be determined by the formula
`i / (2 * n_neigh - i)`, where `i` is a number of nearest neighbours
between 0 and `n_neigh`.

## Examples

``` r
set.seed(2024)
# create an artificial pca embedding
pca_embedding <- matrix(
    c(runif(100 * 10), runif(100 * 10, min = 3, max = 4)),
    nrow = 200, byrow = TRUE
)
rownames(pca_embedding) <- as.character(1:200)
colnames(pca_embedding) <- paste("PC", 1:10)

# calculate the nn adjacency matrix
nn_matrix <- getNNmatrix(
    RANN::nn2(pca_embedding, k = 5)$nn.idx,
    5,
    0,
    -1
)$nn

get_highest_prune_param(nn_matrix, 5)$prune_value
#> [1] 0
```
