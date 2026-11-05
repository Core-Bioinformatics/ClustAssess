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
  aggregate ECC values (for example `element_consistency()` and unordered
  merging).
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
- Leiden with `leidenbase` (`clustering_algorithm = 6`): a worker reuses the
  igraph and edge weights of the last adjacency matrix it clustered for the
  next clusterings of the same matrix, rather than building them for each
  clustering and start. Each worker keeps a copy of that matrix's data and
  its igraph.

Results are those of 1.2.0, except for the fixes under "Correctness and
robustness" and these implementation effects:

- ECS sums and products run in another order, so values can differ in their
  last bits (about 1e-16). A merge threshold equal to a score at that precision
  can then fall on the other side of it.
- Flat ECS calls on character labels run in native code, which creates R's
  random number state (`.Random.seed`) when none exists. They draw no random
  numbers.
- Serial similarity matrices no longer warn that `%dopar%` runs sequentially,
  and their errors no longer start with `task 1 failed`.

## API and options

- Added `make_pca_matrix_processing()` for explicit exact PCA and optional
  fixed-seed `irlba`-backed approximate PCA. The factory uses the input's rank
  limits and ambient BLAS settings, whereas the default processor caps components
  at half the feature count, with a minimum feature-derived cap of one, and
  controls BLAS threads when `RhpcBLASctl` is available. Exact
  factory processing is therefore not necessarily equivalent to the default.
- Added `compute_umap_embeddings` to `assess_feature_stability()` and
  `automatic_stability_assessment()` so callers can skip feature UMAP output
  when only stability scores are required.
- Added `return_nn_idx` to `assess_nn_stability()` and
  `automatic_stability_assessment()` to reduce result size when stored neighbor
  indices are not needed.
- Added opt-in `return_partition_seeds` to `assess_nn_stability()` to expose
  the existing merged representative's seed and RNG kinds for each UMAP graph
  type and neighbor count, without retaining all embeddings or neighbor matrices.

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
- Flat membership vectors containing `NA` or `NaN` now stop with a clear error
  as soon as a function that compares partitions receives them (including
  duplicate merging), instead of failing in native code; the internal
  clustering object is no longer created from them. `merge_partitions()` and
  `element_sim_matrix()` still accept a list of one such partition, as they
  compare nothing. Empty-string labels remain valid; missing observations must
  be handled explicitly before comparison.
- Hardened native batched ECC contingency-table dimension and index arithmetic
  against integer overflow. This does not make intrinsically large tables fit
  available memory.
- Fixed exact duplicate merging for numeric relabellings whose label ranges
  contain unused values.
- Fixed the `avg_agreement` value of `merge_partitions()` (frequency and
  agreement ordering): after duplicate or threshold merging, it now uses the
  frequencies of the merged partitions. Before, it could use the frequencies of
  the input records, which gave wrong values (for example a low value to a
  partition found in most repetitions) and a wrong order of tied partitions.
  Ordering by agreement, and breaking frequency ties by it, can change.
- `merge_partitions()` with frequency or agreement ordering now gives
  `avg_agreement` 1 when every repetition found the same partition. Before,
  that result had no value, or kept a value given with its input. A single
  repetition has no value.
- Fixed exact duplicate merging of factors to ignore unused factor levels.
- Fixed `merge_partitions()` to forward `return_ecs_matrix` and `check_ties`
  through pre-merged and nested partition inputs.
- Fixed `merge_partitions(return_ecs_matrix = TRUE)` so the full returned ECS
  matrix follows partition order after frequency tie-breaking.
- Fixed UMAP argument handling in feature stability, graph stability, and
  connected-component assessment to honor `umap_arguments$n_neighbors` instead
  of silently replacing it with 15. The default remains 15, capped at the cell
  count minus one. Results can change for previously ignored explicit settings.
- Fixed UMAP graph-stability handling for custom seed sequences and for
  repeated neighbor counts, which stopped with an error. A repeated count is
  assessed once, as in PCA graph stability.
- Fixed NN-only and SNN-only UMAP graph-stability output names. PCA NN-only
  graph stability (`graph_type = 0`) no longer returns extra `PCA_snn`
  summaries.
- Fixed Leiden worker dependency selection in `assess_clustering_stability()`
  to use `clustering_algorithm` rather than the auxiliary arguments list.
- Fixed Leiden with `leidenbase` (`clustering_algorithm = 6`):
  - its results are named `Leidenbase`; they were stored under `NA`, which the
    Shiny app and `get_clusters_from_clustassess_object()` could not read;
  - it uses Seurat's partition type, `RBConfigurationVertexPartition`
    (modularity with a resolution parameter, on the scale of Louvain), unless
    `clustering_arguments` sets `partition_type`. `leidenbase`'s own default,
    the Constant Potts Model, reads the resolution on another scale and put
    almost every cell in its own cluster at the usual resolutions;
  - `assess_feature_stability()`, `assess_nn_stability()` and the matching
    steps of `automatic_stability_assessment()` accept it, as
    `assess_clustering_stability()` did.
  - with several starts (`num_starts` above 1 in `clustering_arguments`), it
    gives the same partitions for the same seeds. Each start now gets a seed
    drawn from the clustering's seed; it had none, so `leidenbase` seeded it
    from the clock and the partitions changed from run to run. R's random
    number state after the clustering is as before.
- Fixed `assess_feature_stability()` with `graph_reduction_type = "UMAP"`,
  also as a step of `automatic_stability_assessment()`: each repetition's UMAP
  now starts from the seed of that repetition, as in `assess_nn_stability()`.
  It started from R's random number state at that time unless
  `umap_arguments` had a `seed`, so runs on several workers gave other results
  each time. Results change where `umap_arguments` had no `seed`.
- Fixed the point shapes of `plot_k_resolution_corresp()`,
  `plot_k_n_partitions()` and the Shiny app for five or more clustering
  methods or configurations: the fifth had no shape, or the plot failed. The
  shapes repeat after the five filled ones. The app's method filters show
  every method by default, `Leiden_matrix` and `Leidenbase` included.
- Fixed the Shiny app's plot of the number of partitions per number of
  clusters with its legend at small text sizes, which stopped with "figure
  margins too large": its legend row now holds the margins of its plots.
- Fixed the Shiny app's plots with a legend (the plot of the number of
  partitions per number of clusters, the boxplots of the clustering and
  feature stability, and the UMAPs) at large text or legend sizes or on
  narrow screens, which stopped with "figure margins too large" when the
  legend left the plot no room: the legend text is then drawn smaller, or the
  plot without its legend when even a small one does not fit. A UMAP colour
  bar at a small legend size, under its UMAP or on its own, which stopped in
  the same way, now fits. Where the plot and its legend had 0.1 inch or more,
  they are drawn as before.
- Fixed `automatic_stability_assessment()` to record the actual embedding
  dimension in `stable_config$n_pcs`, including for custom matrix processors.
- Unified feature-count normalization between feature and automatic assessment,
  fixing cache misses for large double-valued counts and fractional counts.
  Counts below one or above the feature-set size, including infinities, use all
  features before integer conversion; missing counts are omitted. Cache
  fallbacks and single-feature processors receive dimension-preserving
  cells-by-features matrices. Default PCA retains one component for a
  one-feature subset.
- Fixed combined PCA/UMAP graph assessment to merge named summaries without
  corrupting neighbor indices, and to pass clustering arguments to both paths.
  A UMAP winner now rebuilds the final graph from its evaluated representative
  seed and worker RNG kinds, including PSOCK RNG streams. Replay restores the
  caller's RNG state. The seed is recorded as `stable_config$umap_seed`, and the
  graph-building embedding is returned when requested. PCA winners keep PCA
  neighbor indices.
- UMAP seed results are now merged in supplied seed order, with or without
  optional provenance. This makes representative selection deterministic with
  respect to worker completion order; tied threshold-merge results can differ
  from previous backend-dependent ordering. The feature-stage UMAPs in
  `feature_stability$embedding_list` stay unchanged; the `umap` of a selected
  configuration can be replaced by the embedding of the winning UMAP graph.
- Made `RhpcBLASctl` and `SharedObject` genuinely optional on default PCA and
  parallel nearest-neighbor paths. Missing packages leave BLAS settings alone
  or use ordinary matrix exports; behavior with the packages installed is retained.
- Exact PCA remains the default. Approximate PCA remains
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
