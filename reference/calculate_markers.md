# Calculate markers

Performs the Wilcoxon rank sum test to identify differentially expressed
genes between two groups of cells.

## Usage

``` r
calculate_markers(
  expression_matrix,
  cells1,
  cells2,
  logfc_threshold = 0,
  min_pct_threshold = 0.1,
  avg_expr_threshold_group1 = 0,
  min_diff_pct_threshold = -Inf,
  rank_matrix = NULL,
  feature_names = NULL,
  used_slot = "data",
  norm_method = "SCT",
  pseudocount_use = 1,
  base = 2,
  adjust_pvals = TRUE,
  check_cells_set_diff = TRUE
)
```

## Arguments

- expression_matrix:

  A matrix of gene expression values having genes in rows and cells in
  columns.

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

- avg_expr_threshold_group1:

  The minimum average expression that a gene should have in the first
  group of cells to be considered as differentially expressed. Defaults
  to `0`.

- min_diff_pct_threshold:

  The minimum difference in the fraction of cells expressing a gene
  between the two cell populations to consider the gene as
  differentially expressed. Defaults to `-Inf`.

- rank_matrix:

  A matrix where the cells are ranked based on their expression levels
  with respect to each gene. Defaults to `NULL`, in which case the
  function will calculate the rank matrix. We recommend calculating the
  rank matrix beforehand and passing it to the function to speed up the
  computation.

- feature_names:

  A vector of gene names. Defaults to `NULL`, in which case the function
  will use the row names of the expression matrix as gene names.

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

- pseudocount_use:

  The pseudocount to add to the expression values when calculating the
  average expression of the genes, to avoid the 0 value for the
  denominator. Defaults to `1`.

- base:

  The base of the logharithm. Defaults to `2`.

- adjust_pvals:

  A logical value indicating whether to adjust the p-values for multiple
  testing using the Bonferonni method. Defaults to `TRUE`.

- check_cells_set_diff:

  A logical value indicating whether to check if thw two cell groups are
  disjoint or not. Defaults to `TRUE`.

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

## Examples

``` r
set.seed(2024)
# create an artificial expression matrix
expr_matrix <- matrix(
    c(runif(100 * 50), runif(100 * 50, min = 3, max = 4)),
    ncol = 200, byrow = FALSE
)
colnames(expr_matrix) <- as.character(1:200)
rownames(expr_matrix) <- paste("feature", 1:50)

calculate_markers(
    expression_matrix = expr_matrix,
    cells1 = 101:200,
    cells2 = 1:100
)
#>                  gene avg_log2FC pct.1 pct.2        p_val    p_val_adj
#> feature 1   feature 1   5.324445     1     1 2.562144e-34 1.281072e-32
#> feature 2   feature 2   5.586295     1     1 2.562144e-34 1.281072e-32
#> feature 3   feature 3   5.273713     1     1 2.562144e-34 1.281072e-32
#> feature 4   feature 4   5.643523     1     1 2.562144e-34 1.281072e-32
#> feature 5   feature 5   5.706965     1     1 2.562144e-34 1.281072e-32
#> feature 6   feature 6   5.487495     1     1 2.562144e-34 1.281072e-32
#> feature 7   feature 7   5.487029     1     1 2.562144e-34 1.281072e-32
#> feature 8   feature 8   5.587873     1     1 2.562144e-34 1.281072e-32
#> feature 9   feature 9   5.578639     1     1 2.562144e-34 1.281072e-32
#> feature 10 feature 10   5.433841     1     1 2.562144e-34 1.281072e-32
#> feature 11 feature 11   5.616236     1     1 2.562144e-34 1.281072e-32
#> feature 12 feature 12   5.467072     1     1 2.562144e-34 1.281072e-32
#> feature 13 feature 13   5.444738     1     1 2.562144e-34 1.281072e-32
#> feature 14 feature 14   5.550933     1     1 2.562144e-34 1.281072e-32
#> feature 15 feature 15   5.464486     1     1 2.562144e-34 1.281072e-32
#> feature 16 feature 16   5.473382     1     1 2.562144e-34 1.281072e-32
#> feature 17 feature 17   5.536340     1     1 2.562144e-34 1.281072e-32
#> feature 18 feature 18   5.553369     1     1 2.562144e-34 1.281072e-32
#> feature 19 feature 19   5.550821     1     1 2.562144e-34 1.281072e-32
#> feature 20 feature 20   5.414216     1     1 2.562144e-34 1.281072e-32
#> feature 21 feature 21   5.339008     1     1 2.562144e-34 1.281072e-32
#> feature 22 feature 22   5.486371     1     1 2.562144e-34 1.281072e-32
#> feature 23 feature 23   5.514801     1     1 2.562144e-34 1.281072e-32
#> feature 24 feature 24   5.626331     1     1 2.562144e-34 1.281072e-32
#> feature 25 feature 25   5.520828     1     1 2.562144e-34 1.281072e-32
#> feature 26 feature 26   5.439109     1     1 2.562144e-34 1.281072e-32
#> feature 27 feature 27   5.510551     1     1 2.562144e-34 1.281072e-32
#> feature 28 feature 28   5.414908     1     1 2.562144e-34 1.281072e-32
#> feature 29 feature 29   5.536102     1     1 2.562144e-34 1.281072e-32
#> feature 30 feature 30   5.471044     1     1 2.562144e-34 1.281072e-32
#> feature 31 feature 31   5.712207     1     1 2.562144e-34 1.281072e-32
#> feature 32 feature 32   5.482001     1     1 2.562144e-34 1.281072e-32
#> feature 33 feature 33   5.396843     1     1 2.562144e-34 1.281072e-32
#> feature 34 feature 34   5.458924     1     1 2.562144e-34 1.281072e-32
#> feature 35 feature 35   5.466319     1     1 2.562144e-34 1.281072e-32
#> feature 36 feature 36   5.521587     1     1 2.562144e-34 1.281072e-32
#> feature 37 feature 37   5.585943     1     1 2.562144e-34 1.281072e-32
#> feature 38 feature 38   5.537862     1     1 2.562144e-34 1.281072e-32
#> feature 39 feature 39   5.531601     1     1 2.562144e-34 1.281072e-32
#> feature 40 feature 40   5.646930     1     1 2.562144e-34 1.281072e-32
#> feature 41 feature 41   5.322801     1     1 2.562144e-34 1.281072e-32
#> feature 42 feature 42   5.709088     1     1 2.562144e-34 1.281072e-32
#> feature 43 feature 43   5.784678     1     1 2.562144e-34 1.281072e-32
#> feature 44 feature 44   5.662688     1     1 2.562144e-34 1.281072e-32
#> feature 45 feature 45   5.541069     1     1 2.562144e-34 1.281072e-32
#> feature 46 feature 46   5.445265     1     1 2.562144e-34 1.281072e-32
#> feature 47 feature 47   5.515937     1     1 2.562144e-34 1.281072e-32
#> feature 48 feature 48   5.621451     1     1 2.562144e-34 1.281072e-32
#> feature 49 feature 49   5.516103     1     1 2.562144e-34 1.281072e-32
#> feature 50 feature 50   5.515016     1     1 2.562144e-34 1.281072e-32
#>            avg_expr_group1
#> feature 1         3.497559
#> feature 2         3.496166
#> feature 3         3.445801
#> feature 4         3.504912
#> feature 5         3.577405
#> feature 6         3.525528
#> feature 7         3.485442
#> feature 8         3.506795
#> feature 9         3.476262
#> feature 10        3.484119
#> feature 11        3.498119
#> feature 12        3.493887
#> feature 13        3.471559
#> feature 14        3.498410
#> feature 15        3.511520
#> feature 16        3.537302
#> feature 17        3.482097
#> feature 18        3.507179
#> feature 19        3.493560
#> feature 20        3.431861
#> feature 21        3.465087
#> feature 22        3.477709
#> feature 23        3.487408
#> feature 24        3.527158
#> feature 25        3.513766
#> feature 26        3.469997
#> feature 27        3.462232
#> feature 28        3.501975
#> feature 29        3.491901
#> feature 30        3.515717
#> feature 31        3.466661
#> feature 32        3.534487
#> feature 33        3.483160
#> feature 34        3.481944
#> feature 35        3.504899
#> feature 36        3.521041
#> feature 37        3.508432
#> feature 38        3.471637
#> feature 39        3.521095
#> feature 40        3.505936
#> feature 41        3.480344
#> feature 42        3.528489
#> feature 43        3.516821
#> feature 44        3.487346
#> feature 45        3.502275
#> feature 46        3.484552
#> feature 47        3.528281
#> feature 48        3.517302
#> feature 49        3.513613
#> feature 50        3.482948
# TODO should be rewritten such that you don't create new matrix objects inside
# just
```
