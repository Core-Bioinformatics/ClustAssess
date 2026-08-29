# Server - Dimensionality reduction module

Creates the backend interface for the dimensionality reduction module
inside the ClustAssess Shiny application.

## Usage

``` r
server_dimensionality_reduction(id, parent_session)
```

## Arguments

- id:

  The id of the module, used to acess the UI elements.

- parent_session:

  The session of the parent module, used to control the tabs of the
  application.

## Note

This function should not be called directly, but in the context of the
app that is created using the `write_shiny_app` function.
