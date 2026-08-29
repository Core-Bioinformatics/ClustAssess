# Server - Graph clustering module

Creates the backend interface for the graph clustering module inside the
ClustAssess Shiny application.

## Usage

``` r
server_graph_clustering(id, feature_choice, parent_session)
```

## Arguments

- id:

  The id of the module, used to acess the UI elements.

- feature_choice:

  A reactive object that contains the chosen configuration from the
  Dimensionality Reduction tab.

- parent_session:

  The session of the parent module, used to control the tabs of the
  application.

## Note

This function should not be called directly, but in the context of the
app that is created using the `write_shiny_app` function.
