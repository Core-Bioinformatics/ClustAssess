# Server - Comparison module

Creates the backend interface for the comparison module inside the
ClustAssess Shiny application.

## Usage

``` r
server_comparisons(id, chosen_config, chosen_method)
```

## Arguments

- id:

  The id of the module, used to acess the UI elements.

- chosen_config:

  A reactive object that contains the chosen configuration from the
  Dimensionality Reduction tab.

- chosen_method:

  A reactive object that contains the chosen method from the Clustering
  tab.

## Note

This function should not be called directly, but in the context of the
app that is created using the `write_shiny_app` function.
