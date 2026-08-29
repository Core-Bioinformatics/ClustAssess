# Create the ClustAssess ShinyApp

Creates the ClustAssess ShinyApp based on the output of the automatic
ClustAssess pipeline. In addition to that, the expression matrix and the
metadata dataframe are provided as input to the ShinyApp. If the
clustassess object is not provided, the function will create the light
version of the ClustAssess ShinyApp, that will not contain the
assessment results. For this case, the `metadata` parameter should
contain two aditional columns named 'UMAP_1" and 'UMAP_2' that will
correspond to the 2D embedding of the cells.

## Usage

``` r
write_shiny_app(
  object,
  metadata = NULL,
  assay_name = NULL,
  clustassess_object,
  project_folder,
  compression_level = 6,
  summary_function = stats::median,
  shiny_app_title = "",
  organism_enrichment = "hsapiens",
  height_ratio = 0.6,
  qualpalr_colorspace = list(h = c(0, 360), s = c(0.1, 0.5), l = c(0.6, 0.85)),
  prompt_feature_choice = TRUE
)

# S3 method for class 'Seurat'
write_shiny_app(
  object,
  metadata = NULL,
  assay_name,
  clustassess_object = NULL,
  project_folder,
  compression_level = 6,
  summary_function = stats::median,
  shiny_app_title = "",
  organism_enrichment = "hsapiens",
  height_ratio = 0.6,
  qualpalr_colorspace = list(h = c(0, 360), s = c(0.1, 0.5), l = c(0.6, 0.85)),
  prompt_feature_choice = TRUE
)

# Default S3 method
write_shiny_app(
  object,
  metadata = NULL,
  assay_name = NULL,
  clustassess_object = NULL,
  project_folder,
  compression_level = 6,
  summary_function = stats::median,
  shiny_app_title = "",
  organism_enrichment = "hsapiens",
  height_ratio = 0.6,
  qualpalr_colorspace = list(h = c(0, 360), s = c(0.1, 0.5), l = c(0.6, 0.85)),
  prompt_feature_choice = TRUE
)
```

## Arguments

- object:

  A Seurat object or an expression matrix

- metadata:

  The metadata dataframe. This parameter will be ignored if the object
  is a Seurat object.

- assay_name:

  The name of the assay to be used to extract the expression matrix from
  the Seurat object. This parameter will be ignored if the object is not
  a Seurat object.

- clustassess_object:

  The output of the ClustAssess automatic pipeline. If the ClustAssess
  object is not provided (NULL), the function will create the light
  version of the ShinyApp, that will not contain the assessment results.

- project_folder:

  The folder where the files will be written

- compression_level:

  The compression level for the h5 files (See \`rhdf5::h5createFile“ for
  more details)

- summary_function:

  The function used for summarizing the stability values; the default is
  `median`

- shiny_app_title:

  The title of the shiny app

- organism_enrichment:

  The organism used for the enrichment analysis; the default is
  `hsapiens`

- height_ratio:

  The ratio of the height of the plot to the height of the browser; the
  default is `0.6`

- qualpalr_colorspace:

  The colorspace used for generating the colors; we use the default as
  defined in the `qualpalr` package

- prompt_feature_choice:

  Should the user be prompted to choose if he wants to continue with the
  selection of features even if it is lower than median sequence depth;
  the default is `TRUE`
