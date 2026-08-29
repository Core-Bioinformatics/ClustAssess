# Server - Landing page module

Creates the backend interface for the landing page module inside the
ClustAssess Shiny application.

## Usage

``` r
server_landing_page(
  id,
  height_ratio,
  dimension,
  parent_session,
  organism = "hsapiens"
)
```

## Arguments

- id:

  The id of the module, used to acess the UI elements.

- height_ratio:

  A reactive object that contains the height ratio of the plots in the
  application (the height of the plot is calculated using the height
  ratio and the height of the webpage).

- dimension:

  A reactive object that contains the dimensions of the webpage.

- parent_session:

  The session of the parent module, used to control the tabs of the
  application.

- organism:

  The organism of the dataset, which will be used in the enrichment
  analysis.

## Note

This function should not be called directly, but in the context of the
app that is created using the `write_shiny_app` function.
