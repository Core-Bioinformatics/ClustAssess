# Changelog

All notable changes to ClustAssess are documented in this file.

This changelog records release contents only.
Earlier release history, including 1.2.0, is retained in [NEWS.md](NEWS.md).

## [1.3.0] - Unreleased

### Algorithm And Performance

- Optimized flat disjoint ECS partition comparison paths used by
  `merge_partitions()` and related consistency scoring.
- Fixed large relabelled duplicate partition merges that could exceed R's
  environment key length limits.
- Avoided unnecessary pairwise ECS matrix construction when callers only need
  aggregate ECC values.
- Reused exact single-cluster-count summary ECC instead of recomputing the same
  weighted consistency result.
- Avoided hidden PSOCK dispatch for non-identity ECS-threshold partition
  merging, keeping that internal path on the native serial calculation.
- Added an adaptive C++ fast path for weighted flat-disjoint ECC calculations,
  retaining the original pairwise path for tiny lists, mixed named/unnamed
  memberships, and unsupported inputs.
- Centralized NN/SNN graph construction from precomputed neighbor indices and
  reduced PCA graph-stability memory pressure by streaming per-neighborhood
  graph matrices.
- Improved final clustering-stability scheduling by distributing
  resolution/seed tasks across workers within each clustering algorithm while
  preserving the existing nested result structure and seed inputs.
- Made forced garbage collection in stability workflows opt-in via
  `options(ClustAssess.force_gc = TRUE)`.

### API And Options

- Added `make_pca_matrix_processing()` for explicit exact PCA processing and
  optional deterministic `irlba`-backed approximate PCA processing.
- Added `compute_umap_embeddings` to `assess_feature_stability()` and
  `automatic_stability_assessment()` so callers can skip feature UMAP output
  when only stability scores are required.
- Added `return_nn_idx` to `assess_nn_stability()` and
  `automatic_stability_assessment()` to reduce result size when stored neighbor
  indices are not needed.
- Frequency-ordered `merge_partitions()` calls no longer attach the
  undocumented `avg_agreement` field unless agreement is needed for ordering or
  tie resolution, or is already available from a requested ECS matrix.

### Correctness And Robustness

- Fixed average ECS for partitions whose numeric label ranges contain unused
  values.
- Fixed exact duplicate merging for numeric relabellings whose label ranges
  contain unused values.
- Fixed `merge_partitions(order_logic = "avg_agreement")` to weight agreement
  values by merged partition frequencies.
- Fixed `merge_partitions()` to forward `return_ecs_matrix` and `check_ties`
  through pre-merged and nested partition inputs.
- Fixed UMAP graph-stability handling for custom seed sequences.
- Fixed NN-only and SNN-only UMAP graph-stability output names.
- Fixed installed-package PSOCK worker execution paths for graph and clustering
  stability workflows.
- Ported the upstream `main` fix selecting the Leiden worker dependency from
  `clustering_algorithm` in `assess_clustering_stability()`, rather than looking
  for an algorithm in the auxiliary `clustering_arguments` list.
- Preserved default PCA and clustering behavior. Approximate PCA remains
  opt-in, uses a fixed seed for deterministic matrix processing, and restores
  the caller RNG state.

### Testing And Benchmarks

- Added correctness and regression coverage for ECS edge cases, PAC scoring,
  graph construction, matrix processing, automatic stability assessment, final
  clustering stability, and PSOCK worker paths.
- Added benchmark and profiling harnesses for ECS, graph construction, matrix
  processing, automatic stability assessment, and clustering-stability worker
  scaling.
- Added a release-branch comparison benchmark harness with optional function
  tracing, `Rprof`/`Rprofmem` output, and saved-result parity checks.
- Added a Linux process-session RSS/PSS sampler so PSOCK memory comparisons
  include the main process and all workers.

### Documentation

- Added portable [performance notes](PERFORMANCE.md) and
  [benchmark instructions](analysis/bench/README.md), including correctness
  boundaries, validation commands, and limitations of the historical evidence.
- Made real-data benchmark paths explicit caller inputs; excluded local agent
  records and raw benchmark artifacts from Git and package builds.
