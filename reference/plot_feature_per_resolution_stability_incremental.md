# Per resolution - Feature Stability Incremental Boxplot

Perform an incremental ECS between two consecutive feature steps. The
ECS values are extracted only for a specified resolution value.

## Usage

``` r
plot_feature_per_resolution_stability_incremental(
  feature_object_list,
  resolution,
  dodge_width = 0.7,
  text_size = 4,
  boxplot_width = 0.4,
  return_df = FALSE
)
```

## Arguments

- feature_object_list:

  An object or a concatenation of objects returned by the
  `assess_feature_stability` method.

- resolution:

  The resolution value for which the ECS will be extracted.

- dodge_width:

  Used for adjusting the horizontal position of the boxplot; the value
  will be passed to the `width` argument of the
  [`ggplot2::position_dodge`](https://ggplot2.tidyverse.org/reference/position_dodge.html)
  method.

- text_size:

  The size of the labels above boxplots.

- boxplot_width:

  Used for adjusting the width of the boxplots; the value will be passed
  to the `width` argument of the
  [`ggplot2::geom_boxplot`](https://ggplot2.tidyverse.org/reference/geom_boxplot.html)
  method.

- return_df:

  If TRUE, the function will return the ECS values as a dataframe.
  Default is `FALSE`.

## Value

A ggplot2 object with ECS distribution will be displayed as a boxplot.
Above each boxplot there will be a pair of numbers representing the two
steps that are compared.

## Examples

``` r
set.seed(2024)
# create an artificial expression matrix
expr_matrix <- matrix(
    c(runif(50 * 10), runif(50 * 10, min = 3, max = 4)),
    nrow = 100, byrow = TRUE
)
rownames(expr_matrix) <- as.character(1:100)
colnames(expr_matrix) <- paste("feature", 1:10)

feature_stability_result <- assess_feature_stability(
    data_matrix = t(expr_matrix),
    feature_set = colnames(expr_matrix),
    steps = c(5, 10),
    feature_type = "feature_name",
    resolution = c(0.1, 0.5),
    n_repetitions = 3,
    umap_arguments = list(
        # the following parameters are used by the umap function
        # and are not mandatory
        n_neighbors = 3,
        approx_pow = TRUE,
        n_epochs = 0,
        init = "random",
        min_dist = 0.3
    ),
    clustering_algorithm = 1
)
plot_feature_per_resolution_stability_incremental(feature_stability_result, 0.1)
```
