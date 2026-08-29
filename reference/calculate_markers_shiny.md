# Calculate markers - Shiny

Performs the Wilcoxon rank sum test to identify differentially expressed
genes between two groups of cells in the shiny context. The method can
be also used outside the shiny context, as long as the expression matrix
is stored in a h5 file.

## Usage

``` r
calculate_markers_shiny(
  cells1,
  cells2,
  logfc_threshold = 0,
  min_pct_threshold = 0.1,
  average_expression_threshold = 0,
  average_expression_group1_threshold = 0,
  min_diff_pct_threshold = -Inf,
  used_slot = "data",
  norm_method = "SCT",
  expression_h5_path = "expression.h5",
  pseudocount_use = 1,
  base = 2,
  verbose = TRUE,
  check_difference = TRUE
)
```

## Arguments

- cells1:

  A vector of cell indices for the first group of cells.

- cells2:

  A vector of cell indices for the second group of cells.

- logfc_threshold:

  The minimum absolute log fold change to consider a gene as
  differentially expressed. Defaults to `0`, meaning all genes are taken
  into considereation.

- min_pct_threshold:

  The minimum fraction of cells expressing a gene form each cell
  population to consider the gene as differentially expressed.
  Increasing the value will speed up the function. Defaults to `0.1`.

- average_expression_threshold:

  The minimum average expression that a gene should have in order to be
  considered as differentially expressed.

- average_expression_group1_threshold:

  The minimum average expression that a gene should have in the first
  group of cells to be considered as differentially expressed. Defaults
  to `0`.

- min_diff_pct_threshold:

  The minimum difference in the fraction of cells expressing a gene
  between the two cell populations to consider the gene as
  differentially expressed. Defaults to `-Inf`.

- used_slot:

  Parameter that provides additional information about the expression
  matrix, whether it was scaled or not. The value of this parameter
  impacts the calculation of the fold change. If `data`, the function
  will calculates the fold change as the fraction between the log value
  of the average of the expression raised to exponential for the two
  cell groups. If `scale.data`, the function will calculate the fold
  change as the fraction between the average of the expression values
  for the two cell groups. Other options will default to calculating the
  fold change as the fraction between the log value of the average of
  the expression values for the two cell groups. Defaults to `data`.

- norm_method:

  The normalization method used to normalize the expression matrix. The
  value of this parameter impacts the calculation of the average
  expression of the genes when `used_slot = "data"`. If `LogNormalize`,
  the log fold change will be calculated as described for the
  `used_slot` parameter. Otherwise, the log fold change will be
  calculated as the fraction between the log value of the average of the
  expression values for the two cell groups. Defaults to `SCT`.

- expression_h5_path:

  The path to the h5 file containing the expression matrix. The h5 file
  should contain the following fields: `expression_matrix`,
  `rank_matrix`, `average_expression`, `genes`. The file path defaults
  to `expression.h5`.

- pseudocount_use:

  The pseudocount to add to the expression values when calculating the
  average expression of the genes, to avoid the 0 value for the
  denominator. Defaults to `1`.

- base:

  The base of the logharithm. Defaults to `2`.

- verbose:

  Whether to print messages about the progress of the function. Defaults
  to TRUE.

- check_difference:

  Whether to perform set difference between the two cells. Defaults to
  TRUE.

## Value

A data frame containing the following columns:

- `gene`: The gene name.

- `avg_log2FC`: The average log fold change between the two cell groups.

- `p_val`: The p-value of the Wilcoxon rank sum test.

- `p_val_adj`: The adjusted p-value of the Wilcoxon rank sum test.

- `pct.1`: The fraction of cells expressing the gene in the first cell
  group.

- `pct.2`: The fraction of cells expressing the gene in the second cell
  group.

- `avg_expr_group1`: The average expression of the gene in the first
  cell group.

- `avg_expr`: The average expression of the gene.
