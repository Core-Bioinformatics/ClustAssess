# Create Clustering Object

Creates a Clustering object from the output of a clustering method.

## Usage

``` r
create_clustering(clustering_result, ...)

# S4 method for class 'numeric'
create_clustering(clustering_result, alpha = 0.9)

# S4 method for class 'integer'
create_clustering(clustering_result, alpha = 0.9)

# S4 method for class 'character'
create_clustering(clustering_result, alpha = 0.9)

# S4 method for class 'factor'
create_clustering(clustering_result, alpha = 0.9)

# S4 method for class 'matrix'
create_clustering(
  clustering_result,
  alpha = 0.9,
  ppr_implementation = "prpack",
  row_normalize = TRUE
)

# S4 method for class 'Matrix'
create_clustering(
  clustering_result,
  alpha = 0.9,
  ppr_implementation = "prpack",
  row_normalize = TRUE
)

# S4 method for class 'hclust'
create_clustering(
  clustering_result,
  alpha = 0.9,
  r = 1,
  rescale_path_type = "max",
  ppr_implementation = "prpack",
  dist_rescaled = FALSE
)
```

## Arguments

- clustering_result:

  The clustering result, either:

  - A numeric/character/factor vector of cluster labels for each
    element.

  - A samples x clusters matrix/Matrix::Matrix of nonzero membership
    values.

  - An hclust object.

- ...:

  This argument is not used.

- alpha:

  A numeric giving the personalized PageRank damping factor; 1 - alpha
  is the restart probability for the PPR random walk.

- ppr_implementation:

  Choose a implementation for personalized page-rank calculation:

  - "prpack": use PPR algorithms in igraph.

  - "power_iteration": use power_iteration method.

- row_normalize:

  Whether to normalize all rows in clustering_result so they sum to one
  before calculating ECS. It is recommended to set this to TRUE, which
  will lead to slightly different ECS values compared to clusim.

- r:

  A numeric hierarchical scaling parameter.

- rescale_path_type:

  A string; rescale the hierarchical height by:

  - "max" : the maximum path from the root.

  - "min" : the minimum path form the root.

  - "linkage" : use the linkage distances in the clustering.

- dist_rescaled:

  A logical: if TRUE, the linkage distances are linearly rescaled to be
  in-between 0 and 1.

## Value

A Clustering object.

## Functions

- `create_clustering(numeric)`: Create Clustering Object from Numeric
  Vector

- `create_clustering(integer)`: Create Clustering Object from Integer
  Vector

- `create_clustering(character)`: Create Clustering Object from
  Character Vector

- `create_clustering(factor)`: Create Clustering Object from Factor
  Vector

- `create_clustering(matrix)`: Create Clustering Object from base matrix

- `create_clustering(Matrix)`: Create Clustering Object from
  Matrix::Matrix

- `create_clustering(hclust)`: Create Clustering Object from hclust
