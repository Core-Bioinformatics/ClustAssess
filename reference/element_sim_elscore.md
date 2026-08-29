# The Element-Centric Clustering Similarity for each Element

Calculates the element-wise element-centric similarity between two
clustering results.

## Usage

``` r
element_sim_elscore(
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

Vector of element-centric similarity between the two clusterings for
each element.

## References

Gates, A. J., Wood, I. B., Hetrick, W. P., & Ahn, Y. Y. (2019).
Element-centric clustering comparison unifies overlaps and hierarchy.
Scientific reports, 9(1), 1-13.
https://doi.org/10.1038/s41598-019-44892-y

## Examples

``` r
km.res <- kmeans(iris[, 1:4], centers = 8)$cluster
hc.res <- hclust(dist(iris[, 1:4]))
element_sim_elscore(km.res, hc.res)
#>          1          2          3          4          5          6          7 
#> 0.35970732 0.28463113 0.27504654 0.28629676 0.34023494 0.14877225 0.24578347 
#>          8          9         10         11         12         13         14 
#> 0.35414515 0.26102546 0.28463113 0.37292373 0.18140078 0.27187915 0.24565854 
#>         15         16         17         18         19         20         21 
#> 0.20410258 0.20410258 0.20482698 0.35970732 0.14877225 0.37418879 0.37070483 
#>         22         23         24         25         26         27         28 
#> 0.37418879 0.20140004 0.36563267 0.18140078 0.27187915 0.36563267 0.35031846 
#>         29         30         31         32         33         34         35 
#> 0.35031846 0.25489307 0.25489307 0.37070483 0.22568160 0.22568160 0.28463113 
#>         36         37         38         39         40         41         42 
#> 0.20413084 0.37632204 0.34023494 0.26102546 0.35414515 0.35993563 0.06290065 
#>         43         44         45         46         47         48         49 
#> 0.26172270 0.36711590 0.36987484 0.28463113 0.37793070 0.28629676 0.37292373 
#>         50         51         52         53         54         55         56 
#> 0.19109146 0.56158716 0.56950919 0.56158716 0.51263990 0.56290102 0.50763058 
#>         57         58         59         60         61         62         63 
#> 0.56950919 0.39088001 0.56290102 0.52643208 0.41153617 0.51868674 0.49906433 
#>         64         65         66         67         68         69         70 
#> 0.55523802 0.51409810 0.57095737 0.50763058 0.53340573 0.54022294 0.51294180 
#>         71         72         73         74         75         76         77 
#> 0.55275532 0.51868674 0.54271092 0.57003772 0.57014260 0.57095737 0.57086373 
#>         78         79         80         81         82         83         84 
#> 0.58021144 0.56221674 0.51409810 0.50674233 0.50674233 0.52168669 0.53530182 
#>         85         86         87         88         89         90         91 
#> 0.50763058 0.57896174 0.57013180 0.54022294 0.52689686 0.51263990 0.50763058 
#>         92         93         94         95         96         97         98 
#> 0.55523802 0.52168669 0.39088001 0.52938883 0.52041953 0.52041953 0.57014260 
#>         99        100        101        102        103        104        105 
#> 0.40000525 0.52938883 0.32101651 0.53373774 0.37552633 0.42059736 0.31830150 
#>        106        107        108        109        110        111        112 
#> 0.38418459 0.43395801 0.37600855 0.31146160 0.38881815 0.42625089 0.55256576 
#>        113        114        115        116        117        118        119 
#> 0.32523715 0.54150766 0.56413558 0.43446073 0.41378213 0.38881815 0.39272284 
#>        120        121        122        123        124        125        126 
#> 0.55356905 0.33853655 0.55087880 0.38418459 0.54279062 0.34649405 0.36750946 
#>        127        128        129        130        131        132        133 
#> 0.54279062 0.54321087 0.31112970 0.36750946 0.37600855 0.38881815 0.31112970 
#>        134        135        136        137        138        139        140 
#> 0.53530182 0.57093378 0.38881815 0.30853205 0.41378213 0.54321087 0.32523715 
#>        141        142        143        144        145        146        147 
#> 0.34166244 0.32523715 0.53373774 0.33853655 0.34166244 0.32523715 0.55256576 
#>        148        149        150 
#> 0.42625089 0.40329848 0.56442420 
```
