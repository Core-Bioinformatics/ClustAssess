## version 1.2.0

---

### Updates
- Add subsetting option for the Shiny heatmaps and bubbleplots. The subset is based on already existing metadata columns.
- Use `leidenbase` package for Leiden clustering instead of `leiden` (aligning with Seurat's changes). Algorithm=6 contains the leidenbase implementation, whereas algorithm=4 and 5 are associated with leiden (for matrix and igraph input objects, respectively).
- Minimise redundant calculations of the adjacency matrix in the automatic assessment.
- Add function that calculates the nn2 index in parallel.

### Fixes
- Ensure the expression matrix or embedding has row names. If not, set them to "cell_1", "cell_2", etc.
- Treat the case when enrichment analysis returns no results.
- Sort the group 2 markers for enrichment analysis based on the absolute logfoldchange.
- Discard the metadata columns with too many unique values for visualisation in the Shiny app.

---

## version 1.1.0

---

### Updates
- Add barplot with Cell count or percentage of metadata in the Shiny context.
- Add option to combine (split) metadata and dynamically create a new metadata column in the Shiny context.
- Add option to calculate the percentage of cells expressing gene above a threshold in the summary table from the Shiny Violin section.
- Add hierarchical plot that shows the relationship between partitions with different number of clusters.
- Add the option to create the ClustAssess app without the need to run the stability assessment (the light version). If the clustassess parameter is NULL, the app will contain only the 'Comparison' tab. In this case, the user should provide the UMAP coordinates in the metadata dataframe.
- Enable live filtering of the marker genes based on all the columns. The filtered table will be used as input for the enrichment analysis.

### Fixes
- Fix the case in `write_object` when the gene variance filtering leaves the chunk with one or zero genes.
- Stop allowing the user to calculate the ECC or perform the merging to a list with less than two partitions by raising an exception.
- Sort the k values numerically in the `merge_resolutions` function.
- Fix the `qualpalr` colour space parameter in the `add_metadata` function.
