# Create monocle object from a ClustAssess shiny app

Use the files generated in the ClustAssess app to create a Monocle
object which has the stable number of clusters.

## Usage

``` r
create_monocle_from_clustassess_app(
  app_folder,
  stable_feature_type,
  stable_feature_set_size,
  stable_clustering_method,
  stable_n_clusters = NULL,
  use_all_genes = FALSE,
  nrows_chunk = 1000,
  refresh_interval = 1
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

- nrows_chunk:

  The size of the chunk (in number of rows) to read from the
  expression.h5 file to form the sparse matrix. Defaults to 1000.

- refresh_interval:

  The periodicity of chunks reading prior to running the garbage
  collector. Defaults to 1 iteration.

## Value

A Monocle object of the expression matrix, having the stable number of
clusters identified by ClustAssess.
