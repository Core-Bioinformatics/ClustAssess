#!/usr/bin/env Rscript

fixture_h5 <- Sys.getenv("CA_ROI_FIXTURE", unset = "")
if (!nzchar(fixture_h5)) {
    stop("Set CA_ROI_FIXTURE to an expression.h5 fixture file.")
}
if (!file.exists(fixture_h5)) {
    stop("Missing fixture: ", fixture_h5)
}

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(x)) y else x
}

args <- commandArgs(FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grepl(file_arg, args)][1])
script_path <- normalizePath(script_path)
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."))

use_installed <- tolower(Sys.getenv("CA_ROI_LOAD_INSTALLED", "false")) %in%
    c("1", "true", "yes")
if (use_installed) {
    suppressPackageStartupMessages(library(ClustAssess))
} else {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
        stop("Package `pkgload` is required unless CA_ROI_LOAD_INSTALLED=true.")
    }
    pkgload::load_all(repo_root, quiet = TRUE)
}

env_int <- function(name, default) {
    value <- Sys.getenv(name, as.character(default))
    as.integer(value)
}

env_num_seq <- function(name, default) {
    value <- Sys.getenv(name, default)
    as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
}

if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Package `rhdf5` is required to run this benchmark.")
}

read_fixture <- function(ngenes, ncells) {
    genes <- rhdf5::h5read(fixture_h5, "genes")
    cells <- rhdf5::h5read(fixture_h5, "cells")
    avg <- rhdf5::h5read(fixture_h5, "average_expression")
    idx <- order(avg, decreasing = TRUE)[seq_len(min(ngenes, length(avg)))]
    cidx <- seq_len(min(ncells, length(cells)))
    expr <- rhdf5::h5read(fixture_h5, "expression_matrix", index = list(idx, cidx))
    rownames(expr) <- genes[idx]
    colnames(expr) <- cells[cidx]
    expr
}

make_inputs <- function(expr) {
    vars <- matrixStats::rowVars(expr)
    list(
        features_sets = list(
            Most_Abundant = rownames(expr),
            Highly_Variable = rownames(expr)[order(vars, decreasing = TRUE)]
        ),
        steps = list(
            Most_Abundant = c(500L, 1000L, min(1500L, nrow(expr))),
            Highly_Variable = c(500L, 1000L, min(1500L, nrow(expr)))
        )
    )
}

register_workers <- function(workers) {
    if (workers <= 1L) {
        foreach::registerDoSEQ()
        return(NULL)
    }
    cl <- parallel::makeCluster(workers)
    doParallel::registerDoParallel(cl)
    cl
}

run_assessment <- function(expr, inputs, workers) {
    cl <- register_workers(workers)
    on.exit({
        if (!is.null(cl)) parallel::stopCluster(cl)
        foreach::registerDoSEQ()
    }, add = TRUE)

    base_mp <- make_pca_matrix_processing(
        actual_npcs = 30L,
        method = "irlba",
        seed = 1L
    )
    mp <- function(dt_mtx, actual_npcs = 30, ...) {
        base_mp(dt_mtx, actual_npcs = actual_npcs, ...)
    }
    formals(mp)$actual_npcs <- 30L

    gc()
    elapsed <- system.time(suppressWarnings(automatic_stability_assessment(
        expression_matrix = expr,
        n_repetitions = env_int("CA_ROI_REPETITIONS", 20L),
        n_neigh_sequence = seq(5L, 50L, 10L),
        resolution_sequence = seq(0.5, 2.5, 0.1),
        features_sets = inputs$features_sets,
        steps = inputs$steps,
        n_top_configs = env_int("CA_ROI_TOP_CONFIGS", 2L),
        matrix_processing = mp,
        umap_arguments = list(
            n_neighbors = 30L,
            n_threads = 1L,
            n_sgd_threads = 1L,
            n_epochs = 10L
        ),
        compute_umap_embeddings = FALSE,
        return_nn_idx = FALSE,
        verbose = FALSE,
        save_temp = FALSE
    )))

    data.frame(
        requested_workers = workers,
        elapsed = unname(elapsed[["elapsed"]]),
        user = unname(elapsed[["user.self"]]),
        sys = unname(elapsed[["sys.self"]]),
        error = NA_character_
    )
}

expr <- read_fixture(
    ngenes = env_int("CA_ROI_GENES", 1500L),
    ncells = env_int("CA_ROI_CELLS", 600L)
)
inputs <- make_inputs(expr)
workers <- as.integer(env_num_seq("CA_ROI_WORKERS", "1"))

results <- do.call(rbind, lapply(workers, function(w) {
    tryCatch(
        run_assessment(expr, inputs, w),
        error = function(e) {
            data.frame(
                requested_workers = w,
                elapsed = NA_real_,
                user = NA_real_,
                sys = NA_real_,
                error = conditionMessage(e)
            )
        }
    )
}))

cat("fixture=", fixture_h5, "\n", sep = "")
cat("slice_genes=", nrow(expr), " slice_cells=", ncol(expr), "\n", sep = "")
cat("settings=repetitions=", env_int("CA_ROI_REPETITIONS", 20L),
    " neigh=5,15,25,35,45 resolutions=21 top_configs=",
    env_int("CA_ROI_TOP_CONFIGS", 2L), " umap=FALSE pca=irlba\n", sep = "")
print(results, row.names = FALSE)
