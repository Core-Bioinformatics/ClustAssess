# Relationship Between Number of Nearest Neighbours and Number of Clusters

Display the distribution of the number of clusters obtained for each
number of neighbours across random seeds.

## Usage

``` r
plot_n_neigh_k_correspondence(nn_object_n_clusters)
```

## Arguments

- nn_object_n_clusters:

  An object or a concatenation of objects returned by the
  `get_nn_importance` method.

## Value

A ggplot2 object with the distributions displayed as boxplots.

## Note

The number of clusters is displayed on a logarithmic scale.

## Examples

``` r
set.seed(2024)
# create an artificial PCA embedding
pca_emb <- matrix(runif(100 * 30), nrow = 100, byrow = TRUE)
rownames(pca_emb) <- as.character(1:100)
colnames(pca_emb) <- paste0("PC_", 1:30)

nn_stability_obj <- assess_nn_stability(
    embedding = pca_emb,
    n_neigh_sequence = c(10, 15, 20),
    n_repetitions = 10,
    graph_reduction_type = "PCA",
    clustering_algorithm = 1
)
plot_n_neigh_k_correspondence(nn_stability_obj)
#> Warning: log-10 transformation introduced infinite values.
```
