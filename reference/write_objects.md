# Write the objects for the ClustAssess ShinyApp

Given the output of the ClustAssess pipeline, the expression matrix and
the metadata, this function creates the files needed for the ClustAssess
ShinyApp. The files are written in the project_folder and are the
following:

- metadata.rds: the metadata file

- stability.h5: contains the stability results

- expression.h5: contains the expression matrix and the rank matrix

## Usage

``` r
write_objects(
  clustassess_object,
  expression_matrix,
  metadata,
  project_folder = ".",
  compression_level = 6,
  chunk_size = 100,
  gene_variance_threshold = 0,
  summary_function = stats::median,
  qualpalr_colorspace = list(h = c(0, 360), s = c(0.1, 0.5), l = c(0.6, 0.85))
)
```

## Arguments

- clustassess_object:

  The output of the ClustAssess automatic pipeline

- expression_matrix:

  The expression matrix

- metadata:

  The metadata

- project_folder:

  The folder where the files will be written

- compression_level:

  The compression level for the h5 files (See \`rhdf5::h5createFile“ for
  more details)

- chunk_size:

  The chunk size for the rank matrix (See `rhdf5::h5createDataset` for
  more details)

- gene_variance_threshold:

  The threshold for the gene variance; genes with variance below this
  threshold will be removed

- summary_function:

  The function used for summarizing the stability values; the default is
  `median`

- qualpalr_colorspace:

  The colorspace used for generating the colors; we use the default as
  defined in the `qualpalr` package

## Value

NULL (the files are written in the project_folder)
