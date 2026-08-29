# Add metadata to ClustAssess ShinyApp

Adds new metadata into the ClustAssess ShinyApp without having to update
the object and re-create the app.

## Usage

``` r
add_metadata(
  app_folder,
  metadata,
  qualpalr_colorspace = list(h = c(0, 360), s = c(0.1, 0.5), l = c(0.6, 0.85))
)
```

## Arguments

- app_folder:

  The folder containing the ClustAssess ShinyApp

- metadata:

  The new metadata to be added. This parameter should be a dataframe
  that follows the same row ordering as the already existing metadata
  from the ClustAssess app.

- qualpalr_colorspace:

  The colorspace to be used for the metadata

## Value

NULL - the metadata object is updated in the app folder
