# Calculate the highest pruning parameter for the SNN graph given Embedding

Given an embedding, the function calculates the highest pruning
parameter for the SNN graph that preserves the connectivity of the
graph.

## Usage

``` r
get_highest_prune_param_embedding(embedding, n_neigh)
```

## Arguments

- embedding:

  A matrix associated with a PCA embedding. Embeddings from other
  dimensionality reduction techniques (such as LSI) can be used.

- n_neigh:

  The number of nearest neighbours.

## Value

The value of the highest pruning parameter.

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

get_highest_prune_param_embedding(pca_embedding, 5)
#> [1] 0
```
