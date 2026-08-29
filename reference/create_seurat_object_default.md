# Create Seurat object

Use a normalized expression matrix and, potentially, an already
generated PCA / UMAP embedding, to create a Seurat object.

## Usage

``` r
create_seurat_object_default(
  normalized_expression_matrix,
  count_matrix = NULL,
  pca_embedding = NULL,
  umap_embedding = NULL,
  metadata_df = NULL,
  verbose = FALSE
)
```

## Arguments

- normalized_expression_matrix:

  The normalized expression matrix having genes on rows and cells on
  columns.

- count_matrix:

  The count matrix having genes on rows and cells on columns. If NULL,
  the normalized_expression_matrix will be used.

- pca_embedding:

  The PCA embedding of the expression matrix. If NULL, the pca will be
  created using the `Seurat` package (default parameters).

- umap_embedding:

  The UMAP embedding of the expression matrix. If NULL, the umap will be
  created using the `Seurat` package (default parameters).

- metadata_df:

  The metadata dataframe having the cell names as rownames. If NULL, a
  dataframe with a single column named `identical_ident` will be
  created.

- verbose:

  Set the level of verbosity of the Seurat functions. Defaults to
  `FALSE`.

## Value

A Seurat object of the expression matrix, having the stable number of
clusters identified by ClustAssess.
