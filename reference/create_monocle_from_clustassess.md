# Create monocle object from a ClustAssess object

Use the object generated using the ClustAssess
`automatic_stability_assessment` function to create a Monocle object
which has the stable number of clusters.

## Usage

``` r
create_monocle_from_clustassess(
  normalized_expression_matrix,
  count_matrix = NULL,
  clustassess_object,
  metadata_df,
  stable_feature_type,
  stable_feature_set_size,
  stable_clustering_method,
  stable_n_clusters = NULL,
  use_all_genes = FALSE
)
```

## Arguments

- normalized_expression_matrix:

  The normalized expression matrix having genes on rows and cells on
  columns.

- count_matrix:

  The count matrix having genes on rows and cells on columns. If NULL,
  the normalized_expression_matrix will be used.

- clustassess_object:

  The output of the `automatic_stability_assessment`.

- metadata_df:

  The metadata dataframe having the cell names as rownames. If NULL, a
  dataframe with a single column named `identical_ident` will be
  created.

- stable_feature_type:

  The feature type which leads to stable clusters.

- stable_feature_set_size:

  The feature size which leads to stable clusters.

- stable_clustering_method:

  The clustering method which leads to stable clusters.

- stable_n_clusters:

  The number of clusters that are stable. If NULL, all the clusters will
  be provided. Defaults to `NULL`.

- use_all_genes:

  A boolean value indicating if the expression matrix should be
  truncated to the genes used in the stability assessment. Defaults to
  `FALSE`.

## Value

A Monocle object of the expression matrix, having the stable number of
clusters identified by ClustAssess.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(2024)
# create an already-transposed artificial expression matrix
expr_matrix <- matrix(
    c(runif(20 * 10), runif(30 * 10, min = 3, max = 4)),
    nrow = 10, byrow = FALSE
)
colnames(expr_matrix) <- as.character(seq_len(ncol(expr_matrix)))
rownames(expr_matrix) <- paste("feature", seq_len(nrow(expr_matrix)))

autom_object <- automatic_stability_assessment(
    expression_matrix = expr_matrix,
    n_repetitions = 3,
    n_neigh_sequence = c(5),
    resolution_sequence = c(0.1, 0.5),
    features_sets = list(
        "set1" = rownames(expr_matrix)
    ),
    steps = list(
        "set1" = c(5, 7)
    ),
    umap_arguments = list(
        # the following parameters have been modified
        # from the default values to ensure that the function
        # will run under 5 seconds
        n_neighbors = 3,
        approx_pow = TRUE,
        n_epochs = 0,
        init = "random",
        min_dist = 0.3
    ),
    n_top_configs = 1,
    algorithms_clustering_assessment = 1,
    save_temp = FALSE,
    verbose = FALSE
)


# uncomment to create the monocle object
# mon_obj <- create_monocle_from_clustassess(
#     normalized_expression_matrix = expr_matrix,
#     clustassess_object = autom_object,
#     metadata = NULL,
#     stable_feature_type = "set1",
#     stable_feature_set_size = "5",
#     stable_clustering_method = "Louvain"
# )
} # }
```
