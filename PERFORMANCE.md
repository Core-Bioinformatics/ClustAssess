# Performance and correctness notes

ClustAssess measures clustering robustness across feature selections, graph
parameters, resolutions, and random seeds. These repeated computations are part
of the assessment, not redundant work to remove indiscriminately. Performance
changes must preserve partitions, frequencies, seed associations, ordering,
names, and consistency scores, except for documented bug fixes and opt-in API
changes.

## Implementation

- Duplicate partitions use compact canonical-membership keys, followed by an
  exact membership check so hash collisions cannot merge different partitions.
  Keys and exact comparisons use the same label equivalence without converting
  numeric labels to text.
- Flat numeric labels are identifiers, not integer cluster indices. Native
  scoring uses a shared compact encoding when needed, preserving distinctions
  between fractional, negative, large, or precision-sensitive values and avoiding
  tables sized from the raw label range. Ordinary pairwise calls incur some
  validation/encoding cost; this correction is not itself a speedup claim.
- `merge_partitions()` computes the ECS matrix and average agreement only when
  needed. Frequency ordering no longer attaches the undocumented
  `avg_agreement` field unless needed for ordering/ties or available from an
  explicitly requested ECS matrix. Nested calls forward the same options.
- Weighted flat-disjoint ECC uses a batched native path from five partitions
  upward. Tiny lists, mixed named/unnamed memberships, and unsupported inputs
  retain the pairwise path. Encoding overhead makes the batch path unsuitable
  for every call size; reference tests cover values and name attributes.
- PCA graph stability constructs one neighborhood's graph matrices at a time.
  SNN-only workers do not retain an unused NN matrix. Shared graph construction
  is checked against the original primitives, including automatic pruning.
- Final clustering distributes resolution/seed jobs within each algorithm and
  explicitly restores their original positions when collecting results. The
  single-cluster-count summary reuses the exact ECC already computed by merging.
- Forced garbage collection is opt-in through
  `options(ClustAssess.force_gc = TRUE)`; normal R garbage collection remains active.

Correctness fixes include large-cell-count ECS arithmetic, numeric-label
equivalence and unused values, frequency-weighted average agreement, full ECS
matrix ordering after ties, nested merge options,
UMAP custom seeds and graph names, and
installed-package PSOCK execution. These are fixes rather than claims of exact
parity with broken baseline behavior; see [NEWS.md](NEWS.md).

UMAP respects an explicit `umap_arguments$n_neighbors`, capped at the cell count
minus one; the default is 15. Earlier versions silently replaced explicit values
with that default. Comparisons using such values can therefore produce different
embeddings or downstream results even with matching seeds.

## Explicit tradeoffs

PCA graph construction, including automatic pruning, runs serially in the main
R process for each neighborhood size. It avoids retaining every neighborhood's
NN/SNN matrices at once, but gives up the previous worker-parallel graph-building
stage. This can trade graph memory for wall-clock time; the net pipeline effect
depends on the workload and worker count.

Final clustering schedules all resolution/seed jobs for an algorithm together.
This provides more jobs to distribute across workers, but retains every returned
membership vector for that algorithm before summarizing resolutions. The previous
schedule retained one resolution's memberships at a time. Peak memory can
therefore increase for many resolutions and seeds. With `verbose = TRUE`, progress
advances during summarization, after all clustering jobs for that algorithm
finish; it is not live worker progress.

Disabling repeated explicit `gc()` calls avoids forced collection pauses, but
can also delay memory reclamation and change peak RSS. Set
`options(ClustAssess.force_gc = TRUE)` when that tradeoff is useful, and measure
both runtime and aggregate memory under the chosen policy.

Exact PCA remains the pipeline default. The explicit
`make_pca_matrix_processing(method = "exact")` factory is not a drop-in replica
of the default processor: for a cells-by-features input it caps components at
`min(requested_pcs, n_cells - 1, n_features)`, while the default processor caps
them at `min(requested_pcs, max(1L, n_features %/% 2))`, retaining one component
for a one-feature input. The factory uses the ambient BLAS
thread settings; the default processor sets BLAS threads to the registered
worker count for PCA and sets them to one afterward.
This thread control applies only when the optional `RhpcBLASctl` package is
installed; otherwise the caller's BLAS settings remain unchanged. Parallel
nearest-neighbor workers use shared matrices when `SharedObject` is available
and ordinary matrix exports otherwise, which can require more worker memory.
Different component counts or numerical results are therefore possible even
with the exact method. Automatic assessment records the actual returned
embedding dimension in `stable_config$n_pcs`.

`make_pca_matrix_processing(method = "irlba")` uses fixed-seed approximate PCA
and restores the caller's RNG state. It also caps components below the feature
count, as required by the truncated SVD. Its embeddings and downstream clusterings
need not match exact PCA.

`compute_umap_embeddings = FALSE` omits stored feature UMAP embeddings. Keep the
default `TRUE` for feature plots, Shiny apps exported with `write_shiny_app()`,
and Monocle conversions with `create_monocle_from_clustassess()` that use these
embeddings. `return_nn_idx = FALSE` omits retained neighbor indices from results.
Both are opt-in output changes; defaults retain those outputs.

With `include_umap_nn_assessment = TRUE`, automatic assessment ranks both PCA
and UMAP graph summaries. If UMAP wins, the existing first merged representative
for that graph type and neighbor count supplies the seed (frequency ordering
with agreement tie-breaking). The pipeline replays only that UMAP with the same
input, recorded worker RNG kinds, and normalized single-thread arguments, then
builds the final graph from its neighbors. This adds one UMAP calculation for a
UMAP winner without retaining
every seed's embedding or neighbor matrix. `stable_config$umap_seed` records the
selected seed; retained `nn_idx` and the configuration's `umap` correspond to that
graph. The replay restores the caller's RNG kinds and state, even when workers
use separate streams such as `clusterSetRNGStream()`.
`compute_umap_embeddings = FALSE` still omits stored embeddings, not the UMAP
calculations needed to assess and build graphs.

The representative-graph choice is a workflow policy, not a uniquely optimal
UMAP embedding. Results are merged in the supplied seed order, without sorting
seed values, including when optional seed metadata is disabled. This makes
tie handling independent of worker completion order. The original feature-stage
plotting embeddings in `feature_stability$embedding_list` remain unchanged;
feature plots use those coordinates with feature-stage scores. Selected
configurations' `umap` coordinates instead serve graph/clustering plots and
Monocle export, and can differ when UMAP wins the graph assessment.

In both entrypoints, feature subset counts below one or above the feature-set
size, including infinities, use the complete feature set before integer
conversion. Remaining counts are truncated, sorted, and deduplicated; missing
values (`NA` and `NaN`) are omitted. This keeps selected feature counts and cached
embeddings consistent for integer, double, and fractional inputs.

Stored feature UMAP embeddings are not used to calculate stability scores.
However, computing them resets and consumes R's global RNG state. Later unseeded
stochastic steps, such as custom `matrix_processing` or UMAP graph reduction,
can therefore differ when toggling `compute_umap_embeddings`. Manage seeds in
those steps if comparable results are required. Exact PCA and processors from
`make_pca_matrix_processing()` do not depend on that intervening RNG state.

No hidden cache is applied to arbitrary `matrix_processing` functions: equal
feature sets need not imply equal ordered inputs, and caller code can have side
effects. Reusing a larger feature set's PCA for smaller subsets changes the
assessment. Further graph caching or result-assembly parallelism needs measured
end-to-end benefit before implementation.

## Measurement requirements and limits

Runtime and memory depend on dataset dimensions, partition diversity, assessment
settings, and worker count. Correctness tests and isolated kernel timings do not
establish end-to-end speedups. Full-scale and multi-dataset performance of this
release remains unvalidated; no general speedup or memory-reduction guarantee is
made.

For a meaningful comparison:

- Force native recompilation and install each revision into a separate library.
  Use matching R, dependency, compiler, and BLAS settings, and ensure workers
  load the same package revision as the main process.
- Keep data, feature ordering, seeds, assessment settings, and requested outputs
  identical. Compare approximate PCA or omitted outputs separately from default
  behavior. Check numerical results as well as timings, allowing only documented
  correctness fixes.
- Measure the complete pipeline before optimizing individual stages. Repeat
  representative cases across dataset sizes and multiple datasets; record exact
  revisions, fixture dimensions, settings, and worker counts with the results.
- Include the main R process and every worker in memory measurements. Report
  aggregate RSS and PSS separately: RSS can count shared pages more than once.
  Allocation totals from `Rprofmem()` are not peak resident memory.

## Validation

Tests cover independent pairwise ECC references, edge-case memberships, graph
primitives, deterministic scheduling, optional-output parity, PCA RNG state,
and installed-package PSOCK paths. Run source tests and a clean package check:

```sh
Rscript -e 'pkgbuild::clean_dll("."); testthat::test_local(".")'
R CMD build --no-manual --no-build-vignettes .
R CMD check --no-manual --no-build-vignettes --install-args="--preclean" ClustAssess_*.tar.gz
```

Use a directory with only the intended tarball for the check command. Source
tests skip installed-package-only worker checks; the tarball check runs them.
Optional Python comparisons require separate environments. Keep object files
during package checks (`--preclean`, without `--clean`) so compiled-code
inspection can run. This differs from clean benchmark installations, which may
remove objects afterward. Benchmark commands and privacy precautions are in
[the benchmark guide](analysis/bench/README.md).

Build the pkgdown website from a clean checkout, as in CI, or a staging tree
containing only intended public files. Its top-level Markdown auto-discovery
does not honor Git ignore rules. Review generated pages and the search index
before publication.
