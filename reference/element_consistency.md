# Element-Wise Consistency Between a Set of Clusterings

Inspect the consistency of a set of clusterings by calculating their
element-wise clustering consistency (also known as element-wise
frustration).

## Usage

``` r
element_consistency(
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

A vector containing the element-wise consistency. If
`calculate_sim_matrix` is set to `TRUE`, the element similarity matrix
will be returned as well.

## References

Gates, A. J., Wood, I. B., Hetrick, W. P., & Ahn, Y. Y. (2019).
Element-centric clustering comparison unifies overlaps and hierarchy.
Scientific reports, 9(1), 1-13.
https://doi.org/10.1038/s41598-019-44892-y

## Examples

``` r
# cluster across 20 random seeds
clustering.list <- lapply(1:20, function(x) kmeans(mtcars, centers = 3)$cluster)
element_consistency(clustering.list)
#>           Mazda RX4       Mazda RX4 Wag          Datsun 710      Hornet 4 Drive 
#>           0.7192531           0.7192531           0.8587655           0.5922761 
#>   Hornet Sportabout             Valiant          Duster 360           Merc 240D 
#>           0.7402256           0.5071650           0.7402256           0.8587655 
#>            Merc 230            Merc 280           Merc 280C          Merc 450SE 
#>           0.8587655           0.7192531           0.7192531           0.6785714 
#>          Merc 450SL         Merc 450SLC  Cadillac Fleetwood Lincoln Continental 
#>           0.6785714           0.6785714           0.6927736           0.6927736 
#>   Chrysler Imperial            Fiat 128         Honda Civic      Toyota Corolla 
#>           0.6927736           0.8587655           0.8587655           0.8587655 
#>       Toyota Corona    Dodge Challenger         AMC Javelin          Camaro Z28 
#>           0.8587655           0.6785714           0.6785714           0.7402256 
#>    Pontiac Firebird           Fiat X1-9       Porsche 914-2        Lotus Europa 
#>           0.6927736           0.8587655           0.8587655           0.8587655 
#>      Ford Pantera L        Ferrari Dino       Maserati Bora          Volvo 142E 
#>           0.7402256           0.7192531           0.7402256           0.8587655 
```
