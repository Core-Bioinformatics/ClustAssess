# PAC Convergence Plot

Plot PAC across iterations for a set of k to assess convergence.

## Usage

``` r
pac_convergence(pac_res, k_plot)
```

## Arguments

- pac_res:

  The data.frame output by consensus_cluster.

- k_plot:

  A vector with values of k to plot.

## Value

A ggplot2 object with the convergence plot. Convergence has been reached
when the lines flatten out across k_plot values. out across

## Examples

``` r
pac.res <- consensus_cluster(iris[, 1:4], k_max = 20)
#> Calculating consensus clustering
pac_convergence(pac.res, k_plot = c(3, 5, 7, 9))
```
