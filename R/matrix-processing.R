#' Create a PCA matrix-processing function
#'
#' Builds a function suitable for the `matrix_processing` argument used by
#' [assess_feature_stability()] and [automatic_stability_assessment()]. The
#' returned function receives a cell-by-feature matrix and returns a
#' cell-by-component embedding.
#'
#' @param actual_npcs Number of principal components to request. The returned
#' function uses this value as its `actual_npcs` default, which can be overridden
#' for each call.
#' @param method PCA implementation. `"exact"` uses [stats::prcomp()].
#' `"irlba"` uses `irlba::prcomp_irlba()` and requires the optional `irlba`
#' package.
#' @param center Passed to the PCA implementation.
#' @param scale. Passed to the PCA implementation.
#' @param seed Random seed used by approximate PCA implementations. The seed is
#' applied inside each approximate PCA call and the caller's RNG state is
#' restored afterward. Exact PCA does not use this seed.
#' @param ... Additional arguments passed to the chosen PCA implementation.
#'
#' @details
#' For a cell-by-feature input, `"exact"` caps the requested component count at
#' `min(actual_npcs, nrow(dt_mtx) - 1, ncol(dt_mtx))`. The `"irlba"` method also
#' caps it at `ncol(dt_mtx) - 1`. The row cap applies even when `center = FALSE`.
#'
#' This factory is an explicit alternative to the default processors in
#' [assess_feature_stability()] and [automatic_stability_assessment()], not a
#' promise of identical output. Those defaults cap requested components at
#' `max(1L, ncol(dt_mtx) %/% 2)` and, when `RhpcBLASctl` is installed,
#' temporarily use the registered foreach worker count for BLAS threads before
#' setting it to one. Without that
#' optional package they leave BLAS settings unchanged. Factory-created processors
#' leave BLAS thread settings unchanged. Choose the same component count and thread settings when
#' comparing implementations; approximate PCA can still produce different results.
#'
#' @return A function with signature `function(dt_mtx, actual_npcs, ...)`.
#' @export
#'
#' @examples
#' matrix_processing <- make_pca_matrix_processing(actual_npcs = 20)
#' \dontrun{
#' matrix_processing <- make_pca_matrix_processing(
#'     actual_npcs = 30,
#'     method = "irlba",
#'     seed = 1
#' )
#' }
make_pca_matrix_processing <- function(actual_npcs = 30,
                                       method = c("exact", "irlba"),
                                       center = TRUE,
                                       scale. = FALSE,
                                       seed = 1L,
                                       ...) {
    if (!is.numeric(actual_npcs) || length(actual_npcs) != 1 || is.na(actual_npcs)) {
        stop("actual_npcs parameter should be a single numeric value")
    }
    actual_npcs <- as.integer(actual_npcs)
    if (actual_npcs < 1) {
        stop("actual_npcs parameter should be at least 1")
    }

    method <- match.arg(method)
    if (!is.numeric(seed) || length(seed) != 1 || is.na(seed)) {
        stop("seed parameter should be a single numeric value")
    }
    seed <- as.integer(seed)
    extra_arguments <- list(...)

    with_preserved_rng_seed <- function(seed, expr) {
        has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        if (has_seed) {
            old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        }

        on.exit({
            if (has_seed) {
                assign(".Random.seed", old_seed, envir = .GlobalEnv)
            } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                rm(".Random.seed", envir = .GlobalEnv)
            }
        }, add = TRUE)

        set.seed(seed)
        force(expr)
    }

    processor <- function(dt_mtx, actual_npcs = 30, ...) {
        if (!is.matrix(dt_mtx) && !inherits(dt_mtx, "Matrix")) {
            stop("dt_mtx parameter should be a matrix")
        }

        npc <- min(as.integer(actual_npcs), nrow(dt_mtx) - 1, ncol(dt_mtx))
        if (method == "irlba") {
            npc <- min(npc, ncol(dt_mtx) - 1)
        }
        if (npc < 1) {
            stop("not enough rows or columns to compute a PCA embedding")
        }

        call_arguments <- c(
            list(x = dt_mtx, center = center, scale. = scale.),
            extra_arguments,
            list(...)
        )

        if (method == "exact") {
            pca <- do.call(stats::prcomp, c(call_arguments, list(rank. = npc)))
        } else {
            if (!requireNamespace("irlba", quietly = TRUE)) {
                stop("Package `irlba` is required when method = 'irlba'")
            }
            pca <- with_preserved_rng_seed(
                seed,
                do.call(irlba::prcomp_irlba, c(call_arguments, list(n = npc)))
            )
        }

        embedding <- pca$x[, seq_len(npc), drop = FALSE]
        rownames(embedding) <- rownames(dt_mtx)
        colnames(embedding) <- paste0("PC_", seq_len(ncol(embedding)))
        embedding
    }

    formals(processor)$actual_npcs <- actual_npcs
    processor
}
