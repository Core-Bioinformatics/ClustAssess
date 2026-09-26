# Benchmarks

Run commands from the repository root with R and the package dependencies
available. Source-tree benchmarks additionally use `pkgload`. Real-data cases
need `rhdf5`; parallel comparisons need `doParallel`; approximate PCA cases
need `irlba`. No dataset or machine-specific environment is bundled.

## Synthetic and optional app fixtures

```sh
Rscript analysis/bench/benchmark_current_behavior.R
```

Optional size controls and a saved result table:

```sh
mkdir -p analysis/bench/results
CA_BENCH_N_CELLS=500 CA_BENCH_N_PARTITIONS=80 CA_BENCH_N_CLUSTERS=12 \
  CA_BENCH_INCLUDE_PIPELINE=true \
  CA_BENCH_OUTPUT=analysis/bench/results/current_behavior.csv \
  Rscript analysis/bench/benchmark_current_behavior.R
```

Supply an app directory explicitly to include a real-data case:

```sh
CA_BENCH_APP_DIR="path/to/app" CA_BENCH_APP_CELLS=500 CA_BENCH_APP_GENES=1000 \
  Rscript analysis/bench/benchmark_current_behavior.R
```

Replace placeholder paths with your own inputs. An app's `expression.h5` contains
`expression_matrix` (features by cells), `genes`, `cells`, and
`average_expression`. The general benchmark can also read `stability.h5` when
present. App fixtures are optional coverage, not unit-test inputs. When `irlba`
is available, an additional PCA/NN/SNN case compares exact and approximate PCA;
these methods are not expected to produce identical embeddings.

The synthetic representative automatic-pipeline fixture is opt-in:

```sh
CA_BENCH_INCLUDE_REPRESENTATIVE_PIPELINE=true \
  CA_BENCH_REP_GENES=400 CA_BENCH_REP_CELLS=200 CA_BENCH_REP_NPCS=20 \
  Rscript analysis/bench/benchmark_current_behavior.R
```

Other controls are `CA_BENCH_REP_STEP_SMALL`, `CA_BENCH_REP_STEP_LARGE`, and
`CA_BENCH_REP_REPETITIONS`. This fixture omits optional feature UMAP embeddings
and retained NN indices to focus on stability assessment and matrix processing.

## App-data automatic-pipeline benchmark

`benchmark_canard_fixture.R` accepts any compatible ClustAssess app fixture and
uses `CA_CANARD_*` options. There is no default dataset location: an unset or
empty `CA_CANARD_APP_DIR` is an error.

```sh
CA_CANARD_APP_DIR="path/to/app" \
  Rscript analysis/bench/benchmark_canard_fixture.R
```

This reads a bounded slice of `expression.h5`, derives abundance- and
variance-ranked feature sets, and profiles reduced automatic assessments.
Defaults are:

- `CA_CANARD_GENES=1500`, `CA_CANARD_CELLS=600`, `CA_CANARD_NPCS=30`
- `CA_CANARD_STEPS=500,1000,1500`, `CA_CANARD_REPETITIONS=3`
- `CA_CANARD_N_NEIGH=5,15,25`, `CA_CANARD_RESOLUTION=0.5,1.5,2.5`
- `CA_CANARD_UMAP_EPOCHS=40`
- `CA_CANARD_INCLUDE_EXACT=true`, `CA_CANARD_INCLUDE_IRLBA=true`,
  `CA_CANARD_INCLUDE_UMAP=true`

`CA_CANARD_PROFILE_DIR` selects the profile directory; `CA_CANARD_OUTPUT` can
save the timing table. Choose a runtime and memory budget before increasing
fixture dimensions or repetition counts.

For targeted post-PCA profiling, use `investigate_post_pca_targets.R` with the
required `CA_ROI_FIXTURE="path/to/expression.h5"`. It supports
`CA_ROI_LOAD_INSTALLED=true` to profile an installed package instead of loading
the source tree. It requires `irlba` and uses `matrixStats` for feature variance.

## Release-branch comparison

Install each source tree into a separate empty library with a forced native
rebuild. `R CMD INSTALL` can otherwise reuse stale `src/*.o` files, invalidating
the comparison. Use matching R/dependency versions and compiler/BLAS settings.

```sh
mkdir -p analysis/bench/libraries/baseline analysis/bench/libraries/current
R CMD INSTALL --preclean --clean --library=analysis/bench/libraries/baseline \
  path/to/baseline-checkout
R CMD INSTALL --preclean --clean --library=analysis/bench/libraries/current .
```

With `rhdf5` and `matrixStats` available, run the same fixture/settings in each
library. `CA_COMPARE_FIXTURE` is required; it has no machine-specific default.

```sh
R_LIBS=analysis/bench/libraries/baseline \
  CA_COMPARE_LIB=analysis/bench/libraries/baseline \
  CA_COMPARE_LABEL=baseline CA_COMPARE_MODE=strict \
  CA_COMPARE_FIXTURE="path/to/expression.h5" \
  Rscript analysis/bench/compare_release_branches.R
```

Repeat with the current library and label, leaving data, API options, and seeds
unchanged. Default `strict` mode uses the common public API. `current_fast` uses
1.3.0 opt-in approximate PCA and output omissions, so it is not an exact-output
comparison with baseline defaults.

Comparison controls:

- `CA_COMPARE_GENES`, `CA_COMPARE_CELLS`: bounded real-data slice size.
- `CA_COMPARE_REPETITIONS`, `CA_COMPARE_N_NEIGH`, `CA_COMPARE_RESOLUTION`,
  `CA_COMPARE_TOP_CONFIGS`, `CA_COMPARE_UMAP_EPOCHS`: assessment settings.
- `CA_COMPARE_WORKERS`: number of PSOCK workers (default: 8). Set this explicitly
  to fit the available memory and CPU budget.
- `CA_COMPARE_TRACE=true`, `CA_COMPARE_TRACE_FUNCTIONS`: function tracing.
- `CA_COMPARE_RPROF=true`, `CA_COMPARE_PROFILE_DIR`: CPU/allocation profiles.
- `CA_COMPARE_SAVE_RESULT`: optional RDS output, with its MD5 recorded in the
  timing table, for separate parity checks; use a path under
  `analysis/bench/results/` after creating that directory. The harness does not
  compare saved results. Load the baseline and current RDS files and compare
  their contents separately; a hash alone is not a numerical parity test.

Use `measure_process_tree.py` on Linux to measure the complete main/worker
process session. Main-process RSS alone is misleading when work moves between
the parent process and workers. The sampler reports peak aggregate RSS/PSS:

```sh
R_LIBS=analysis/bench/libraries/current \
  CA_COMPARE_LIB=analysis/bench/libraries/current \
  CA_COMPARE_LABEL=current CA_COMPARE_FIXTURE="path/to/expression.h5" \
  python3 analysis/bench/measure_process_tree.py -- \
  Rscript analysis/bench/compare_release_branches.R
```

For output comparisons, distinguish documented bug fixes and optional
`avg_agreement` metadata from numerical regressions. See
[the performance notes](../../PERFORMANCE.md) for those boundaries and the
measurement requirements. Repeat important cases before making release
performance claims; single runs are useful for scouting, not broad guarantees.

## Weighted ECC and profiling

`benchmark_weighted_ecc.R` compares installed dispatch, an independent pairwise
reference, encoding alone, and the native kernel on pre-encoded input. Controls
include `CA_ECC_LIB`, `CA_ECC_INPUT`, `CA_ECC_ELEMENTS`, `CA_ECC_CLUSTERINGS`,
`CA_ECC_CLUSTERS`, `CA_ECC_SIM_MATRIX`, `CA_ECC_ITERATIONS`, and
`CA_ECC_CALLS_PER_ITERATION`. Without `CA_ECC_INPUT` it generates synthetic data.
Parity checks include legacy names and mixed-name fallback behavior.

General benchmarks use base `system.time()`, one explicit pre-case `gc()`, and
`gcFirst = FALSE`. Set `CA_BENCH_PROFILE=true` for CPU/allocation diagnostics.
Profiles default to `analysis/bench/profiles/`; configure their location,
sampling interval, and allocation threshold with `CA_BENCH_PROFILE_DIR`,
`CA_BENCH_PROFILE_INTERVAL`, and `CA_BENCH_PROFILE_MEMORY_THRESHOLD`.

## Local-data boundary

Keep fixture data in `analysis/bench/fixtures/`, benchmark libraries in
`analysis/bench/libraries/`, and raw results/profiles in `results/` or `profiles/`
under this directory. Those directories are excluded from Git and package
tarballs. The scripts do not upload data. Custom output locations outside these
directories are not automatically protected.

Raw logs, profiles, saved objects, and process-sampler output can contain paths,
commands, or data identifiers. Review and sanitize them before deliberately
publishing results.
