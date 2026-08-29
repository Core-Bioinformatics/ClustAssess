# Server - Sandbox module

Creates the backend interface for the sandbox module inside the
ClustAssess Shiny application.

## Usage

``` r
server_sandbox(id)
```

## Arguments

- id:

  The id of the module, used to acess the UI elements.

## Note

This function should not be called directly, but in the context of the
app that is created using the `write_shiny_app` function.
