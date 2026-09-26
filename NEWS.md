# ClustAssess 1.3.0 (unreleased)

## Algorithm and performance

These are implementation changes, not a measured end-to-end speedup or
memory-reduction guarantee. See `PERFORMANCE.md` in the source repository for
tradeoffs and measurement requirements.

- Reworked flat disjoint ECS partition comparison paths used by
  `merge_partitions()` and related consistency scoring.
- Added compact label-invariant keys for duplicate partition matching, with
  the same label equivalence as the exact comparison and no numeric-to-text
  rounding. Exact comparisons guard against hash collisions.
- Avoided unnecessary pairwise ECS matrix construction when callers only need
  aggregate ECC values.
- Reused exact single-cluster-count summary ECC instead of recomputing the same
  weighted consistency result.
- Avoided hidden PSOCK dispatch for non-identity ECS-threshold partition
  merging, keeping that internal path on the native serial calculation.
- Added an adaptive C++ fast path for weighted flat-disjoint ECC calculations,
  retaining the original pairwise path for tiny lists, mixed named/unnamed
  memberships, and unsupported inputs.
- Centralized NN/SNN graph construction from precomputed neighbor indices. PCA
  graph stability now builds one neighborhood's graph matrices at a time in the
  main process rather than building all neighborhoods in parallel on workers.
- Changed final clustering-stability scheduling to distribute resolution/seed
  tasks across workers within each algorithm, preserving nested results and
  seed inputs. All memberships for that algorithm are retained before
  summarization; progress advances after those jobs finish.
- Made forced garbage collection in stability workflows opt-in via
  `options(ClustAssess.force_gc = TRUE)`. Avoiding forced pauses can also delay
  memory reclamation; normal R garbage collection remains active.

## API and options

- Added `make_pca_matrix_processing()` for explicit exact PCA and optional
  fixed-seed `irlba`-backed approximate PCA. The factory uses the input's rank
  limits and ambient BLAS settings, whereas the unchanged default processor
  caps components at half the feature count and controls BLAS threads. Exact
  factory processing is therefore not necessarily equivalent to the default.
- Added `compute_umap_embeddings` to `assess_feature_stability()` and
  `automatic_stability_assessment()` so callers can skip feature UMAP output
  when only stability scores are required.
- Added `return_nn_idx` to `assess_nn_stability()` and
  `automatic_stability_assessment()` to reduce result size when stored neighbor
  indices are not needed.
- Frequency-ordered `merge_partitions()` calls no longer attach the
  undocumented `avg_agreement` field unless agreement is needed for ordering or
  tie resolution, or is already available from a requested ECS matrix.

## Correctness and robustness

- Fixed integer overflow in `element_sim()` when the cell count times a cluster
  size exceeds the integer range, which could produce negative, infinite, or
  inflated similarity values. The optimized similarity-matrix and threshold
  merging paths also use this corrected calculation.
- Fixed average ECS for partitions whose numeric label ranges contain unused
  values.
- Made flat numeric labels consistent across scalar/per-element ECS, reference
  comparisons with `element_agreement()`, pairwise and batched ECC, similarity
  matrices, and partition merging. Fractional,
  negative, large, and precision-sensitive labels retain their equality pattern
  instead of being truncated or rounded through integer/text conversion.
  Internal compact encoding also avoids allocating tables from raw label ranges.
- Flat membership vectors containing `NA` or `NaN` now produce a clear error,
  including singleton and duplicate merge shortcuts. Empty-string labels remain
  valid; missing observations must be handled explicitly before comparison.
- Hardened native batched ECC contingency-table dimension and index arithmetic
  against integer overflow. This does not make intrinsically large tables fit
  available memory.
- Fixed exact duplicate merging for numeric relabellings whose label ranges
  contain unused values.
- Fixed `merge_partitions(order_logic = "avg_agreement")` to weight agreement
  values by merged partition frequencies.
- Fixed `merge_partitions()` to forward `return_ecs_matrix` and `check_ties`
  through pre-merged and nested partition inputs.
- Fixed `merge_partitions(return_ecs_matrix = TRUE)` so the full returned ECS
  matrix follows partition order after frequency tie-breaking.
- Fixed UMAP argument handling in feature stability, graph stability, and
  connected-component assessment to honor `umap_arguments$n_neighbors` instead
  of silently replacing it with 15. The default remains 15, capped at the cell
  count minus one. Results can change for previously ignored explicit settings.
- Fixed UMAP graph-stability handling for custom seed sequences.
- Fixed NN-only and SNN-only UMAP graph-stability output names.
- Fixed installed-package PSOCK worker execution paths for graph and clustering
  stability workflows.
- Fixed Leiden worker dependency selection in `assess_clustering_stability()`
  to use `clustering_algorithm` rather than the auxiliary arguments list.
- Fixed `automatic_stability_assessment()` to record the actual embedding
  dimension in `stable_config$n_pcs`, including for custom matrix processors.
- Preserved default PCA and clustering behavior. Approximate PCA remains
  opt-in, uses a fixed seed for deterministic matrix processing, and restores
  the caller RNG state.

## Testing and benchmarks

- Added correctness and regression coverage for ECS edge cases, PAC scoring,
  graph construction, matrix processing, automatic stability assessment, final
  clustering stability, and PSOCK worker paths.
- Added benchmark and profiling harnesses for ECS, graph construction, matrix
  processing, automatic stability assessment, and clustering-stability worker
  scaling.
- Added a release-branch comparison benchmark harness with optional function
  tracing, `Rprof`/`Rprofmem` output, and saved results (with MD5) for separate
  parity checks. The harness does not itself compare saved results.
- Added a Linux process-session RSS/PSS sampler so PSOCK memory comparisons
  include the main process and all workers.
- Require testthat 3.2.0 or newer for regression-test binding mocks, and use a
  supported Matrix coercion in clustering test fixtures.

## Documentation

- Added portable `PERFORMANCE.md` and `analysis/bench/README.md` guides with
  correctness boundaries, runtime/memory tradeoffs, validation commands, and
  requirements for comparable measurements.
- Documented downstream UMAP consumers and RNG considerations for stochastic
  steps when optional feature UMAP computation is disabled.
- Updated the pipeline and dimensionality-reduction vignettes for explicit
  UMAP neighbor counts, processor contracts, retained outputs, and actual
  embedding dimensions. Separate dependency/data setup from rendering, clean up
  example worker processes, and fix composite QC figure construction.
- Made real-data benchmark paths explicit caller inputs; excluded raw benchmark
  artifacts from Git and package builds.

# ClustAssess 1.2.0

---

## Updates

- Add subsetting option for the Shiny heatmaps and bubbleplots. The subset is based on already existing metadata columns.
- Use `leidenbase` package for Leiden clustering instead of `leiden` (aligning with Seurat's changes). Algorithm=6 contains the leidenbase implementation, whereas algorithm=4 and 5 are associated with leiden (for matrix and igraph input objects, respectively).
- Minimise redundant calculations of the adjacency matrix in the automatic assessment.
- Add function that calculates the nn2 index in parallel.
- Add an ordering of the points in the metadata shiny plots based on the order of the input subgroups.
- Scale the colour values of gene expression in ShinyApp based on provided cap value.
- Control the number of values shown on the continuous UMAP legend.
- Add option to download multiple genes in separate UMAP files.
- Removes `ggtext` dependency from the package (used in the shiny gene UMAP plots).

## Fixes

- Ensure the expression matrix or embedding has row names. If not, set them to "cell_1", "cell_2", etc.
- Treat the case when enrichment analysis returns no results.
- Sort the group 2 markers for enrichment analysis based on the absolute logfoldchange.
- Discard the metadata columns with too many unique values for visualisation in the Shiny app.
- When creating a Seurat object, the number of PCs is set to the minimum between 30 and half the number of genes. This is to avoid errors when the number of genes is low.

---

# ClustAssess 1.1.0

---

## Updates

- Add barplot with Cell count or percentage of metadata in the Shiny context.
- Add option to combine (split) metadata and dynamically create a new metadata column in the Shiny context.
- Add option to calculate the percentage of cells expressing gene above a threshold in the summary table from the Shiny Violin section.
- Add hierarchical plot that shows the relationship between partitions with different number of clusters.
- Add the option to create the ClustAssess app without the need to run the stability assessment (the light version). If the clustassess parameter is NULL, the app will contain only the 'Comparison' tab. In this case, the user should provide the UMAP coordinates in the metadata dataframe.
- Enable live filtering of the marker genes based on all the columns. The filtered table will be used as input for the enrichment analysis.

## Fixes

- Fix the case in `write_object` when the gene variance filtering leaves the chunk with one or zero genes.
- Stop allowing the user to calculate the ECC or perform the merging to a list with less than two partitions by raising an exception.
- Sort the k values numerically in the `merge_resolutions` function.
- Fix the `qualpalr` colour space parameter in the `add_metadata` function.
