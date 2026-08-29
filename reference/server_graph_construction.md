# Server - Graph construction module

Creates the backend interface for the graph construction module inside
the ClustAssess Shiny application.

## Usage

``` r
server_graph_construction(id, chosen_config)
```

## Arguments

- id:

  The id of the module, used to acess the UI elements.

- chosen_config:

  A reactive object that contains the chosen configuration from the
  Dimensionality Reduction tab.

## Note

This function should not be called directly, but in the context of the
app that is created using the `write_shiny_app` function.
