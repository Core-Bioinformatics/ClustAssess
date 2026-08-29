# Merge Partitions from different Resolutions

Merge partitions obtained with different resolution values. The
partitions will be grouped based on the number of clusters. The
identical partitions will be merged into a single partition by updating
the frequency using the `merge_partitions` method.

## Usage

``` r
merge_resolutions(res_obj)
```

## Arguments

- res_obj:

  A list associated to a configuration field from the object returned by
  the `assess_clustering_importance` method.

## Value

A list having one field assigned to each number of clusters. A number of
cluster will contain a list of all merged partitions. To avoid
duplicates, `merged_partitions` with threshold 1 is applied.
