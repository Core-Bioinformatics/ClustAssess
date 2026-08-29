# Assess the stability for Graph Building Parameters

Evaluates clustering stability when changing the values of different
parameters involved in the graph building step, namely the base
embedding, the graph type and the number of neighbours.

## Usage

``` r
assess_nn_stability(
  embedding,
  n_neigh_sequence,
  n_repetitions = 100,
  seed_sequence = NULL,
  graph_reduction_type = "PCA",
  ecs_thresh = 1,
  graph_type = 2,
  prune_value = -1,
  clustering_algorithm = 1,
  clustering_arguments = list(),
  umap_arguments = list()
)
```

## Arguments

- embedding:

  A matrix associated with a PCA embedding. Embeddings from other
  dimensionality reduction techniques (such as LSI) can be used.

- n_neigh_sequence:

  A sequence of the number of nearest neighbours.

- n_repetitions:

  The number of repetitions of applying the pipeline with different
  seeds; ignored if seed_sequence is provided by the user.

- seed_sequence:

  A custom seed sequence; if the value is NULL, the sequence will be
  built starting from 1 with a step of 100.

- graph_reduction_type:

  The graph reduction type, denoting if the graph should be built on
  either the PCA or the UMAP embedding.

- ecs_thresh:

  The ECS threshold used for merging similar clusterings.

- graph_type:

  Argument indicating whether the graph should be unweighted (0),
  weighted (1) or both (2).

- prune_value:

  Argument indicating whether to prune the SNN graph. If the value is 0,
  the graph won't be pruned. If the value is between 0 and 1, the edges
  with weight under the pruning value will be removed. If the value is
  -1, the highest pruning value will be calculated automatically and
  used.

- clustering_algorithm:

  An index indicating which community detection algorithm will be used:
  Louvain (1), Louvain refined (2), SLM (3) or Leiden (4). More details
  can be found in the Seurat's `FindClusters` function.

- clustering_arguments:

  A list of arguments that will be passed to the clustering algorithm.
  See the `FindClusters` function in Seurat for more details.

- umap_arguments:

  Additional arguments passed to the the
  [`uwot::umap`](https://jlmelville.github.io/uwot/reference/umap.html)
  method.

## Value

A list having three fields:

- `n_neigh_k_corresp` - list containing the number of the clusters
  obtained by running the pipeline multiple times with different seed,
  number of neighbours and graph type (weighted vs unweigted)

- `n_neigh_ec_consistency` - list containing the EC consistency of the
  partitions obtained at multiple runs when changing the number of
  neighbours or the graph type

- `n_different_partitions` - the number of different partitions obtained
  by each number of neighbours

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
plot_n_neigh_ecs(nn_stability_obj)
```
