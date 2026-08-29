# Create Seurat object from a ClustAssess shiny app

Use the files generated in the ClustAssess app to create a Seurat object
which has the stable number of clusters.

## Usage

``` r
create_seurat_object_from_clustassess_app(
  app_folder,
  stable_feature_type,
  stable_feature_set_size,
  stable_clustering_method,
  stable_n_clusters = NULL,
  use_all_genes = FALSE
)
```

## Arguments

- app_folder:

  Path pointing to the folder containing a ClustAssess app.

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

A Seurat object of the expression matrix, having the stable number of
clusters identified by ClustAssess.
