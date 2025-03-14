# Unreleased

### Updates
- Add barplot with Cell count or percentage of metadata in the Shiny context.
- Add option to combine (split) metadata and dynamically create a new metadata column in the Shiny context.
- Add option to calculate the percentage of cells expressing gene above a threshold in the summary table from the Shiny Violin section.

### Fixes
- Fix the case in `write_object` when the gene variance filtering leaves the chunk with one or zero genes.