# The Clustering Class

A class containing relevant data for comparing clusterings, including
the affinity matrix for the Clustering.

## Slots

- `names`:

  A character vector of element names; will be 1:n_elements if no names
  were available when creating the Clustering object.

- `n_elements`:

  A numeric giving the number of elements.

- `is_hierarchical`:

  A logical indicating whether the clustering is hierarchical or flat.

- `is_disjoint`:

  A logical indicating whether the clustering is disjoint or
  overlapping.

- `alpha`:

  A numeric giving the personalized PageRank damping factor; 1 - alpha
  is the restart probability for the PPR random walk.

- `r`:

  A numeric hierarchical scaling parameter.

- `elm2clu_dict`:

  A list giving the clusters each element is a member of.

- `clu2elm_dict`:

  A list giving the element members of each cluster.

- `affinity_matrix`:

  A Matrix containing the personalized pagerank equilibrium
  distribution.
