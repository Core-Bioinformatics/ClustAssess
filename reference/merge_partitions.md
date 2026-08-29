# Merge Partitions

Merge flat disjoint clusterings whose pairwise ECS score is above a
given threshold. The merging is done using a complete linkage approach.

## Usage

``` r
merge_partitions(
  partition_list,
  ecs_thresh = 1,
  order_logic = c("freq", "avg_agreement", "none"),
  return_ecs_matrix = FALSE,
  check_ties = TRUE
)
```

## Arguments

- partition_list:

  A list of flat disjoint membership vectors.

- ecs_thresh:

  A numeric: the ecs threshold.

- order_logic:

  Variable indicating the method of ordering the partitions. It can take
  these three values:

  - "freq": order the partitions based on their frequencies. The
    partition with the highest frequency will be the first on the list
    (default).

  - "avg_agreement": order the partitions based on their average
    agreement index. The average agreement index of a partition is
    calculated as the mean of the ECS scores between that partition and
    the other partitions from the list. The partition with the highest
    agreement will be the first on the list.

  - "none": do not perform any ordering (not recommended). If selected,
    the average agreement scores will not be calculated.

- return_ecs_matrix:

  A logical: if TRUE, the function will add the ECS matrix to the return
  list. Defaults to FALSE.

- check_ties:

  A logical value that indicates whether to check for ties in the
  highest frequency partitions or not. If TRUE, the function will put at
  the first position the partition that has the highest similarity with
  the other partitions. Defaults to `FALSE`.

## Value

a list of the merged partitions, together with their associated ECC
score. If `return_ecs_matrix` is set to TRUE, the function will also
return the ECS matrix.

## Examples

``` r
initial_list <- list(c(1, 1, 2), c(2, 2, 2), c("B", "B", "A"))
merge_partitions(initial_list, 1)
#> $partitions
#> $partitions[[1]]
#> $partitions[[1]]$mb
#> [1] 1 1 2
#> 
#> $partitions[[1]]$freq
#> [1] 2
#> 
#> $partitions[[1]]$avg_agreement
#> [1] 0.2777778
#> 
#> 
#> $partitions[[2]]
#> $partitions[[2]]$mb
#> [1] 2 2 2
#> 
#> $partitions[[2]]$freq
#> [1] 1
#> 
#> $partitions[[2]]$avg_agreement
#> [1] 0.2777778
#> 
#> 
#> 
#> $ecc
#> [1] 0.7777778 0.7777778 0.5555556
#> 
```
