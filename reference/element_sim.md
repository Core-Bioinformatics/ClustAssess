# The Element-Centric Clustering Similarity

Calculates the average element-centric similarity between two clustering
results

## Usage

``` r
element_sim(
  clustering1,
  clustering2,
  alpha = 0.9,
  r_cl1 = 1,
  rescale_path_type_cl1 = "max",
  ppr_implementation_cl1 = "prpack",
  dist_rescaled_cl1 = FALSE,
  row_normalize_cl1 = TRUE,
  r_cl2 = 1,
  rescale_path_type_cl2 = "max",
  ppr_implementation_cl2 = "prpack",
  dist_rescaled_cl2 = FALSE,
  row_normalize_cl2 = TRUE
)
```

## Arguments

- clustering1:

  The first clustering result, which can be one of:

  - A numeric/character/factor vector of cluster labels for each
    element.

  - A samples x clusters matrix/Matrix::Matrix of nonzero membership
    values.

  - An hclust object.

- clustering2:

  The second clustering result, which can be one of:

  - A numeric/character/factor vector of cluster labels for each
    element.

  - A samples x clusters matrix/Matrix::Matrix of nonzero membership
    values.

  - An hclust object.

- alpha:

  A numeric giving the personalized PageRank damping factor; 1 - alpha
  is the restart probability for the PPR random walk.

- r_cl1:

  A numeric hierarchical scaling parameter for the first clustering.

- rescale_path_type_cl1:

  A string; rescale the hierarchical height of the first clustering by:

  - "max" : the maximum path from the root.

  - "min" : the minimum path form the root.

  - "linkage" : use the linkage distances in the clustering.

- ppr_implementation_cl1:

  Choose a implementation for personalized page-rank calculation for the
  first clustering:

  - "prpack": use PPR algorithms in igraph.

  - "power_iteration": use power_iteration method.

- dist_rescaled_cl1:

  A logical: if TRUE, the linkage distances of the first clustering are
  linearly rescaled to be in-between 0 and 1.

- row_normalize_cl1:

  Whether to normalize all rows in the first clustering so they sum to
  one before calculating ECS. It is recommended to set this to TRUE,
  which will lead to slightly different ECS values compared to clusim.

- r_cl2:

  A numeric hierarchical scaling parameter for the second clustering.

- rescale_path_type_cl2:

  A string; rescale the hierarchical height of the second clustering by:

  - "max" : the maximum path from the root.

  - "min" : the minimum path form the root.

  - "linkage" : use the linkage distances in the clustering.

- ppr_implementation_cl2:

  Choose a implementation for personalized page-rank calculation for the
  second clustering:

  - "prpack": use PPR algorithms in igraph.

  - "power_iteration": use power_iteration method.

- dist_rescaled_cl2:

  A logical: if TRUE, the linkage distances of the second clustering are
  linearly rescaled to be in-between 0 and 1.

- row_normalize_cl2:

  Whether to normalize all rows in the second clustering so they sum to
  one before calculating ECS. It is recommended to set this to TRUE,
  which will lead to slightly different ECS values compared to clusim.

## Value

The average element-wise similarity between the two Clusterings.

## Examples

``` r
km.res <- kmeans(mtcars, centers = 3)$cluster
hc.res <- hclust(dist(mtcars))
element_sim(km.res, hc.res)
#> [1] 0.4444444
```
