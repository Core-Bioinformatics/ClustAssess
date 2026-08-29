# Element-Wise Average Agreement Between a Set of Clusterings

Inspect how consistently of a set of clusterings agree with a reference
clustering by calculating their element-wise average agreement.

## Usage

``` r
element_agreement(
  reference_clustering,
  clustering_list,
  alpha = 0.9,
  r = 1,
  rescale_path_type = "max",
  ppr_implementation = "prpack",
  dist_rescaled = FALSE,
  row_normalize = TRUE
)
```

## Arguments

- reference_clustering:

  The reference clustering, that each clustering in clustering_list is
  compared to. It can be either:

  - A numeric/character/factor vector of cluster labels for each
    element.

  - A samples x clusters matrix/Matrix::Matrix of nonzero membership
    values.

  - An hclust object.

- clustering_list:

  The list of clustering results, each of which is either:

  - A numeric/character/factor vector of cluster labels for each
    element.

  - A samples x clusters matrix/Matrix::Matrix of nonzero membership
    values.

  - An hclust object.

- alpha:

  A numeric giving the personalized PageRank damping factor; 1 - alpha
  is the restart probability for the PPR random walk.

- r:

  A numeric hierarchical scaling parameter.

- rescale_path_type:

  A string; rescale the hierarchical height by:

  - "max" : the maximum path from the root.

  - "min" : the minimum path form the root.

  - "linkage" : use the linkage distances in the clustering.

- ppr_implementation:

  Choose a implementation for personalized page-rank calculation:

  - "prpack": use PPR algorithms in igraph.

  - "power_iteration": use power_iteration method.

- dist_rescaled:

  A logical: if TRUE, the linkage distances are linearly rescaled to be
  in-between 0 and 1.

- row_normalize:

  Whether to normalize all rows in clustering_result so they sum to one
  before calculating ECS. It is recommended to set this to TRUE, which
  will lead to slightly different ECS values compared to clusim.

## Value

A vector containing the element-wise average agreement.

## References

Gates, A. J., Wood, I. B., Hetrick, W. P., & Ahn, Y. Y. (2019).
Element-centric clustering comparison unifies overlaps and hierarchy.
Scientific reports, 9(1), 1-13.
https://doi.org/10.1038/s41598-019-44892-y

## Examples

``` r
# perform k-means clustering across 20 random seeds
reference.clustering <- iris$Species
clustering.list <- lapply(1:20, function(x) kmeans(iris[, 1:4], centers = 3)$cluster)
element_agreement(reference.clustering, clustering.list)
#>   [1] 0.8980000 0.8020000 0.8020000 0.8020000 0.8980000 0.8980000 0.8980000
#>   [8] 0.8980000 0.8020000 0.8020000 0.8980000 0.8980000 0.8020000 0.8020000
#>  [15] 0.8980000 0.8980000 0.8980000 0.8980000 0.8980000 0.8980000 0.8980000
#>  [22] 0.8980000 0.8980000 0.8980000 0.8020000 0.8020000 0.8980000 0.8980000
#>  [29] 0.8980000 0.8020000 0.8020000 0.8980000 0.8980000 0.8980000 0.8020000
#>  [36] 0.8980000 0.8980000 0.8980000 0.8020000 0.8980000 0.8980000 0.8020000
#>  [43] 0.8020000 0.8980000 0.8980000 0.8020000 0.8980000 0.8020000 0.8980000
#>  [50] 0.8980000 0.6856855 0.6856855 0.1717500 0.6856855 0.6856855 0.6856855
#>  [57] 0.6856855 0.5659355 0.6856855 0.6856855 0.5659355 0.6856855 0.6856855
#>  [64] 0.6856855 0.6856855 0.6856855 0.6856855 0.6856855 0.6856855 0.6856855
#>  [71] 0.6856855 0.6856855 0.6856855 0.6856855 0.6856855 0.6856855 0.6856855
#>  [78] 0.1717500 0.6856855 0.6856855 0.6856855 0.6856855 0.6856855 0.6856855
#>  [85] 0.6856855 0.6856855 0.6856855 0.6856855 0.6856855 0.6856855 0.6856855
#>  [92] 0.6856855 0.6856855 0.5659355 0.6856855 0.6856855 0.6856855 0.6856855
#>  [99] 0.5659355 0.6856855 0.6602500 0.3143145 0.6602500 0.6602500 0.6602500
#> [106] 0.6602500 0.3143145 0.6602500 0.6602500 0.6602500 0.6602500 0.6602500
#> [113] 0.6602500 0.3143145 0.3143145 0.6602500 0.6602500 0.6602500 0.6602500
#> [120] 0.3143145 0.6602500 0.3143145 0.6602500 0.3143145 0.6602500 0.6602500
#> [127] 0.3143145 0.3143145 0.6602500 0.6602500 0.6602500 0.6602500 0.6602500
#> [134] 0.3143145 0.6602500 0.6602500 0.6602500 0.6602500 0.3143145 0.6602500
#> [141] 0.6602500 0.6602500 0.3143145 0.6602500 0.6602500 0.6602500 0.3143145
#> [148] 0.6602500 0.6602500 0.3143145
```
