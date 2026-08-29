# Clustering Method Stability Facet Plot

Display the distribution of the EC consistency for each clustering
method and each resolution value on a given embedding The `all` field of
the object returned by the `get_clustering_difference_object` method is
used.

## Usage

``` r
plot_clustering_difference_facet(
  clust_object,
  embedding,
  low_limit = 0,
  high_limit = 1,
  grid = TRUE
)
```

## Arguments

- clust_object:

  An object returned by the `assess_clustering_stability` method.

- embedding:

  An embedding (only the first two dimensions will be used for
  visualization).

- low_limit:

  The lowest value of ECC that will be displayed on the embedding.

- high_limit:

  The highest value of ECC that will be displayed on the embedding.

- grid:

  Boolean value indicating whether the facet should be a grid (where
  each row is associated with a resolution value and each column with a
  clustering method) or a wrap.

## Value

A ggplot2 object. \#TODO should export

## Examples

``` r
# FIXME fix the examples
# set.seed(2021)
# # create an artificial PCA embedding
# pca_embedding <- matrix(runif(100 * 30), nrow = 100)
# rownames(pca_embedding) <- as.character(1:100)
# colnames(pca_embedding) <- paste0("PCA_", 1:30)

# adj_matrix <- Seurat::FindNeighbors(pca_embedding,
#     k.param = 10,
#     nn.method = "rann",
#     verbose = FALSE,
#     compute.SNN = FALSE
# )$nn
# clust_diff_obj <- assess_clustering_stability(
#     graph_adjacency_matrix = adj_matrix,
#     resolution = c(0.5, 1),
#     n_repetitions = 10,
#     algorithm = 1:2,
#     verbose = FALSE
# )
# plot_clustering_difference_facet(clust_diff_obj, pca_embedding)
```
