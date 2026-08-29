# Get the vector of colours from a palette

Using the `paletteer` package, this function retrieves a vector of
colours from a specified palette. The function will look for both
discrete and continuous palettes. If the palette is not found, a default
option will be used.

## Usage

``` r
get_colour_vector_from_palette(
  palette_name,
  is_inverse = FALSE,
  placeholder = "viridis::viridis"
)
```

## Arguments

- palette_name:

  The name of the palette to retrieve. The naming should follow the
  guidelines described in the `paletteer` package.

- is_inverse:

  Logical. If `TRUE`, the colours will be reversed.

- placeholder:

  The default palette to use if the specified palette is not found. The
  default is "viridis::viridis".

## Value

A vector of colours from the specified palette. If the palette is not
found, a default palette will be used. If `paletter` is not installed,
an error will be thrown.
