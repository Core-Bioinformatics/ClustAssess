# Choose stable clusters based on ECC and frequency

Filter the list of clusters obtained by the automatic ClustAssess
pipeline using the ECC and frequency thresholds. The ECC threshold is
meant to filter out the partitions that are highly sensitive to the
change of the random seed, while the purpose of the frequency threshold
is to assure a statistical significance of the inferred stability.

## Usage

``` r
choose_stable_clusters(
  clusters_list,
  ecc_threshold = 0.9,
  freq_threshold = 30,
  summary_function = mean
)
```

## Arguments

- clusters_list:

  List of clusters obtained from the
  `get_clusters_from_clustassess_object` function.

- ecc_threshold:

  Minimum ECC value to consider a cluster as stable. Default is 0.9.

- freq_threshold:

  Minimum total frequency of the partitions to consider. Default is 30.

- summary_function:

  Function to summarize the ECC values. Default is `mean`. To match the
  results from the ClustAssess Shiny App, use `median`.

## Value

A list of stable clusters that satisfy the ECC and frequency.
