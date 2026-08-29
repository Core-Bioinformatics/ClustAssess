# Extract config-specific clusters from a ClustAssess object

Given the output of the `automatic_stability_assessment` function,
extract the clusters that are specific to a particular configuration of
feature type, feature size, clustering method and, optionally, the
number of clusters.

## Usage

``` r
get_clusters_from_clustassess_object(
  clustassess_object,
  feature_type = NULL,
  feature_size = NULL,
  clustering_method = NULL,
  nclusters = NULL
)
```

## Arguments

- clustassess_object:

  Output of the `automatic_stability_assessment`.

- feature_type:

  Type of feature used for dimensionality reduction. If `NULL`, it will
  select the first available feature.

- feature_size:

  Size of the feature set used for clustering. If `NULL`, it will select
  the first available feature size.

- clustering_method:

  Clustering method used. If `NULL`, it will select the first available
  clustering method.

- nclusters:

  Number of clusters to extract. If `NULL`, all clusters are returned.

## Value

A list of clusters that are specific to the given configuration. Each
number of cluster will contain the list of partitions with that specific
k and the ECC value indicating the overall stability of k.
