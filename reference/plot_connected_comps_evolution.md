# Relationship Between Number of Nearest Neighbours and Graph Connectivity

Display the distribution of the number connected components obtained for
each number of neighbours across random seeds.

## Usage

``` r
plot_connected_comps_evolution(nn_conn_comps_object)
```

## Arguments

- nn_conn_comps_object:

  An object or a concatenation of objects returned by the
  `get_nn_conn_comps` method.

## Value

A ggplot2 object with boxplots for the connected component
distributions.

## Note

The number of connected components is displayed on a logarithmic scale.

## Examples

``` r
set.seed(2024)
# create an artificial PCA embedding
pca_emb <- matrix(runif(100 * 30), nrow = 100, byrow = TRUE)
rownames(pca_emb) <- as.character(1:100)
colnames(pca_emb) <- paste0("PCA_", 1:30)

nn_conn_comps_obj <- get_nn_conn_comps(
    embedding = pca_emb,
    n_neigh_sequence = c(2, 5),
    n_repetitions = 3,
    # arguments that are passed to the uwot function
    umap_arguments = list(
        min_dist = 0.3,
        metric = "cosine"
    )
)
plot_connected_comps_evolution(nn_conn_comps_obj)
#> Warning: log-10 transformation introduced infinite values.
```
