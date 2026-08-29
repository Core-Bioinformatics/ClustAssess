# Create monocle object

Use a normalized expression matrix and, potentially, an already
generated PCA / UMAP embedding, to create a Monocle object.

## Usage

``` r
create_monocle_default(
  normalized_expression_matrix,
  count_matrix = NULL,
  pca_embedding = NULL,
  umap_embedding = NULL,
  metadata_df = NULL
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
  created using the `monocle3` package (default parameters).

- umap_embedding:

  The UMAP embedding of the expression matrix. If NULL, the umap will be
  created using the `monocle3` package (default parameters).

- metadata_df:

  The metadata dataframe having the cell names as rownames. If NULL, a
  dataframe with a single column named `identical_ident` will be
  created.

## Value

A Monocle object of the expression matrix, having the stable number of
clusters identified by ClustAssess.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(2024)
# create an already-transposed artificial expression matrix
expr_matrix <- matrix(
    c(runif(20 * 10), runif(30 * 10, min = 3, max = 4)),
    nrow = 10, byrow = FALSE
)
colnames(expr_matrix) <- as.character(seq_len(ncol(expr_matrix)))
rownames(expr_matrix) <- paste("feature", seq_len(nrow(expr_matrix)))

# uncomment to create the monocle object
mon_obj <- create_monocle_default(
    normalized_expression_matrix = expr_matrix,
    pca_emb = NULL,
    umap_emb = NULL,
    metadata_df = NULL
)
} # }
```
