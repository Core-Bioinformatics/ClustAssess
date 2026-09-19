# Performance and correctness notes

ClustAssess measures clustering robustness across feature selections, graph
parameters, resolutions, and random seeds. These repeated computations are part
of the assessment, not redundant work to remove indiscriminately. Performance
changes must preserve partitions, frequencies, seed associations, ordering,
names, and consistency scores, except for documented bug fixes and opt-in API
changes.

## Included changes

- Duplicate partitions use compact canonical-membership keys, followed by an
  exact membership check so hash collisions cannot merge different partitions.
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

Correctness fixes include numeric labels with unused values, frequency-weighted
average agreement, nested merge options, UMAP custom seeds and graph names, and
installed-package PSOCK execution. These are fixes rather than claims of exact
parity with broken baseline behavior; see [CHANGELOG.md](CHANGELOG.md).

## Explicit tradeoffs

Exact PCA remains the default. `make_pca_matrix_processing(method = "irlba")`
uses fixed-seed approximate PCA and restores the caller's RNG state, but its
embeddings and downstream clusterings need not match exact PCA.

`compute_umap_embeddings = FALSE` omits feature plotting embeddings, and
`return_nn_idx = FALSE` omits retained neighbor indices from results. Both are
opt-in output changes; defaults retain those outputs.

No hidden cache is applied to arbitrary `matrix_processing` functions: equal
feature sets need not imply equal ordered inputs, and caller code can have side
effects. Reusing a larger feature set's PCA for smaller subsets changes the
assessment. Further graph caching or result-assembly parallelism needs measured
end-to-end benefit before implementation.

## Evidence and limits

Historical July 2026 comparisons used forced clean native builds with matching
compiler settings, the same data/settings/seeds, and memory across the main R
process plus eight PSOCK workers. On one private 3,000-feature by 2,000-cell
slice, the then-current optimization working tree took 79.785s versus 95.067s
for baseline `efba325`, with about 15% lower aggregate peak RSS. The fixture is
not distributed, and these single-run observations are not general speedup
guarantees or a reproducible public-dataset result.

Earlier apparent regressions involved a reused slow native object and
inconsistent main-process versus worker-memory accounting. This is why clean
builds and complete process-memory measurement are required. July timings were
not repeated after integration with baseline `5744cf0`; full-scale and
multi-dataset performance remain unvalidated.

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
