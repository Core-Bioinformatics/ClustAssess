## version 1.2.0

---

### Updates
- Add subsetting option for the Shiny heatmaps and bubbleplots. The subset is based on already existing metadata columns.

### Fixes

---

## version 1.1.0

---

### Updates
- Add barplot with Cell count or percentage of metadata in the Shiny context.
- Add option to combine (split) metadata and dynamically create a new metadata column in the Shiny context.
- Add option to calculate the percentage of cells expressing gene above a threshold in the summary table from the Shiny Violin section.
- Add hierarchical plot that shows the relationship between partitions with different number of clusters.
- Add the option to create the ClustAssess app without the need to run the stability assessment (the light version). If the clustassess parameter is NULL, the app will contain only the 'Comparison' tab. In this case, the user should provide the UMAP coordinates in the metadata dataframe.

### Fixes
- Fix the case in `write_object` when the gene variance filtering leaves the chunk with one or zero genes.
- Stop allowing the user to calculate the ECC or perform the merging to a list with less than two partitions by raising an exception.
- Sort the k values numerically in the `merge_resolutions` function.
