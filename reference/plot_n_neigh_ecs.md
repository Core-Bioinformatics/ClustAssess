# Graph construction parameters - ECC facet

Display, for all configurations consisting in different number of
neighbours, graph types and base embeddings, the EC Consistency of the
partitions obtained over multiple runs on an UMAP embedding.

## Usage

``` r
plot_n_neigh_ecs(nn_ecs_object, boxplot_width = 0.5)
```

## Arguments

- nn_ecs_object:

  An object or a concatenation of objects returned by the
  `get_nn_importance` method.

- boxplot_width:

  Used for adjusting the width of the boxplots; the value will be passed
  to the `width` argument of the
  [`ggplot2::geom_boxplot`](https://ggplot2.tidyverse.org/reference/geom_boxplot.html)
  method.

## Value

A ggplot2 object.

## Examples

``` r
set.seed(2024)
# create an artificial PCA embedding
pca_emb <- matrix(runif(100 * 10), nrow = 100, byrow = TRUE)
rownames(pca_emb) <- as.character(1:100)
colnames(pca_emb) <- paste0("PC_", 1:10)

nn_stability_obj <- assess_nn_stability(
    embedding = pca_emb,
    n_neigh_sequence = c(10, 15, 20),
    n_repetitions = 5,
    graph_reduction_type = "PCA",
    clustering_algorithm = 1
)
plot_n_neigh_ecs(nn_stability_obj)
```
