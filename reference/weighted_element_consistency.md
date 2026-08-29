# Weighted Element-Centric Consistency

Calculate the weighted element-centric consistency of a set of
clusterings. The weights are used to give more importance to some
clusterings over others.

## Usage

``` r
weighted_element_consistency(
  clustering_list,
  weights = NULL,
  calculate_sim_matrix = FALSE
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

- weights:

  A numeric vector of weights for each clustering in `clustering_list`.
  If `NULL`, then all weights will be equal to 1. Defaults to `NULL`.

- calculate_sim_matrix:

  A logical value that indicates whether to calculate the similarity
  matrix or not along with the consistency score. Defaults to `FALSE`.

## Value

A vector containing the weighted element-wise consistency. If
`calculate_sim_matrix` is set to `TRUE`, the element similarity matrix
will be returned as well.

## Note

The weighted ECC will be calculated as \\\displaystyle \frac{\sum\_{i}
\sum\_{j} w_i w_j ECS(i, j)}{\sum\_{i} w_i}\\

## Examples

``` r
# cluster across 20 random seeds
clustering_list <- lapply(1:20, function(x) kmeans(mtcars, centers = 3)$cluster)
weights <- sample(1:10, 20, replace = TRUE)
weighted_element_consistency(clustering_list, weights = weights)
#>           Mazda RX4       Mazda RX4 Wag          Datsun 710      Hornet 4 Drive 
#>           0.6750269           0.6750269           0.8437837           0.5971647 
#>   Hornet Sportabout             Valiant          Duster 360           Merc 240D 
#>           0.7748918           0.5612399           0.7748918           0.8437837 
#>            Merc 230            Merc 280           Merc 280C          Merc 450SE 
#>           0.8437837           0.6750269           0.6750269           0.6750297 
#>          Merc 450SL         Merc 450SLC  Cadillac Fleetwood Lincoln Continental 
#>           0.6750297           0.6750297           0.7391709           0.7391709 
#>   Chrysler Imperial            Fiat 128         Honda Civic      Toyota Corolla 
#>           0.7391709           0.8437837           0.8437837           0.8437837 
#>       Toyota Corona    Dodge Challenger         AMC Javelin          Camaro Z28 
#>           0.8437837           0.6750297           0.6750297           0.7748918 
#>    Pontiac Firebird           Fiat X1-9       Porsche 914-2        Lotus Europa 
#>           0.7391709           0.8437837           0.8437837           0.8437837 
#>      Ford Pantera L        Ferrari Dino       Maserati Bora          Volvo 142E 
#>           0.7748918           0.6750269           0.7748918           0.8437837 
```
