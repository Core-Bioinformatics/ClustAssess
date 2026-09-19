#!/usr/bin/env Rscript

app_dir <- Sys.getenv("CA_CANARD_APP_DIR", unset = "")
if (!nzchar(app_dir)) {
    stop("Set CA_CANARD_APP_DIR to an app directory containing expression.h5.")
}
expr_h5 <- file.path(app_dir, "expression.h5")
if (!file.exists(expr_h5)) {
    stop("Missing expression.h5 in CA_CANARD_APP_DIR: ", app_dir)
}

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grepl(file_arg, args)][1])
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."))

if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Package `pkgload` is required to run this benchmark.")
}
if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Package `rhdf5` is required to run this benchmark.")
}

pkgload::load_all(repo_root, quiet = TRUE)

env_flag <- function(name, default = FALSE) {
    value <- Sys.getenv(name, if (default) "true" else "false")
    tolower(value) %in% c("1", "true", "yes")
}

env_int <- function(name, default) {
    as.integer(Sys.getenv(name, as.character(default)))
}

safe_filename <- function(label) {
    gsub("[^A-Za-z0-9_.-]+", "_", label)
}

row_variance <- function(x) {
    if (requireNamespace("matrixStats", quietly = TRUE)) {
        return(matrixStats::rowVars(x))
    }
    means <- rowMeans(x)
    rowSums((x - means)^2) / pmax(ncol(x) - 1, 1)
}

summarize_rprofmem <- function(memory_file) {
    lines <- readLines(memory_file, warn = FALSE)
    allocation_lines <- grepl("^[0-9]+[[:space:]]*:", lines)
    allocation_bytes <- suppressWarnings(as.numeric(
        sub("^([0-9]+)[[:space:]]*:.*", "\\1", lines[allocation_lines])
    ))

    data.frame(
        allocation_count = length(allocation_bytes),
        allocated_bytes = sum(allocation_bytes, na.rm = TRUE),
        stringsAsFactors = FALSE
    )
}

run_case <- function(label, expr, profile_dir, profile_interval = 0.01) {
    gc()
    dir.create(profile_dir, recursive = TRUE, showWarnings = FALSE)
    profile_file <- file.path(profile_dir, paste0(safe_filename(label), ".Rprof"))
    memory_file <- file.path(profile_dir, paste0(safe_filename(label), ".Rprofmem"))

    utils::Rprof(profile_file, interval = profile_interval, memory.profiling = TRUE)
    utils::Rprofmem(memory_file, threshold = 0)
    elapsed <- system.time(force(expr), gcFirst = FALSE)[["elapsed"]]
    utils::Rprofmem(NULL)
    utils::Rprof(NULL)

    profile_total_file <- file.path(profile_dir, paste0(safe_filename(label), "_by_total.csv"))
    profile_self_file <- file.path(profile_dir, paste0(safe_filename(label), "_by_self.csv"))
    profile_summary <- tryCatch(
        utils::summaryRprof(profile_file, memory = "both"),
        error = function(e) NULL
    )
    if (!is.null(profile_summary)) {
        if (!is.null(profile_summary$by.total) && nrow(profile_summary$by.total) > 0) {
            utils::write.csv(profile_summary$by.total, profile_total_file)
        }
        if (!is.null(profile_summary$by.self) && nrow(profile_summary$by.self) > 0) {
            utils::write.csv(profile_summary$by.self, profile_self_file)
        }
    }

    memory_summary <- summarize_rprofmem(memory_file)
    data.frame(
        case = label,
        elapsed_seconds = elapsed,
        profile_file = profile_file,
        profile_total_file = profile_total_file,
        profile_self_file = profile_self_file,
        profile_memory_file = memory_file,
        allocation_count = memory_summary$allocation_count,
        allocated_bytes = memory_summary$allocated_bytes,
        stringsAsFactors = FALSE
    )
}

profile_dir <- Sys.getenv(
    "CA_CANARD_PROFILE_DIR",
    file.path(repo_root, "analysis", "bench", "profiles", "canard_fixture")
)
profile_interval <- as.numeric(Sys.getenv("CA_CANARD_PROFILE_INTERVAL", "0.01"))

n_genes <- env_int("CA_CANARD_GENES", 1500)
n_cells <- env_int("CA_CANARD_CELLS", 600)
n_pcs <- env_int("CA_CANARD_NPCS", 30)
n_repetitions <- env_int("CA_CANARD_REPETITIONS", 3)
n_epochs <- env_int("CA_CANARD_UMAP_EPOCHS", 40)
steps <- as.integer(strsplit(Sys.getenv("CA_CANARD_STEPS", "500,1000,1500"), ",")[[1]])
n_neigh_sequence <- as.integer(strsplit(Sys.getenv("CA_CANARD_N_NEIGH", "5,15,25"), ",")[[1]])
resolution_sequence <- as.numeric(strsplit(Sys.getenv("CA_CANARD_RESOLUTION", "0.5,1.5,2.5"), ",")[[1]])

genes <- rhdf5::h5read(expr_h5, "genes")
cells <- rhdf5::h5read(expr_h5, "cells")
average_expression <- rhdf5::h5read(expr_h5, "average_expression")
n_genes <- min(n_genes, length(genes))
n_cells <- min(n_cells, length(cells))
selected_gene_idx <- order(average_expression, decreasing = TRUE)[seq_len(n_genes)]
selected_cell_idx <- seq_len(n_cells)

expression_matrix <- rhdf5::h5read(
    expr_h5,
    "expression_matrix",
    index = list(selected_gene_idx, selected_cell_idx)
)
rownames(expression_matrix) <- genes[selected_gene_idx]
colnames(expression_matrix) <- cells[selected_cell_idx]

valid_steps <- steps[steps >= 2 & steps <= nrow(expression_matrix)]
if (length(valid_steps) == 0) {
    stop("No requested CA_CANARD_STEPS fit within the selected gene count.")
}

most_abundant <- rownames(expression_matrix)
highly_variable <- rownames(expression_matrix)[order(row_variance(expression_matrix), decreasing = TRUE)]
feature_sets <- list(
    Most_Abundant = most_abundant,
    Highly_Variable = highly_variable
)
step_sets <- list(
    Most_Abundant = valid_steps,
    Highly_Variable = valid_steps
)

common_arguments <- list(
    expression_matrix = expression_matrix,
    n_repetitions = n_repetitions,
    n_neigh_sequence = n_neigh_sequence,
    resolution_sequence = resolution_sequence,
    features_sets = feature_sets,
    steps = step_sets,
    n_top_configs = 2L,
    umap_arguments = list(
        n_neighbors = 30L,
        min_dist = 0.3,
        metric = "cosine",
        n_threads = 1L,
        n_sgd_threads = 1L,
        n_epochs = n_epochs
    ),
    return_nn_idx = FALSE,
    verbose = FALSE,
    save_temp = FALSE
)

cases <- list()
if (env_flag("CA_CANARD_INCLUDE_EXACT", TRUE)) {
    exact_processor <- make_pca_matrix_processing(actual_npcs = n_pcs, method = "exact")
    cases[["canard_fixture_exact_no_umap"]] <- list(
        matrix_processing = exact_processor,
        compute_umap_embeddings = FALSE
    )
}

if (env_flag("CA_CANARD_INCLUDE_IRLBA", TRUE)) {
    irlba_processor <- make_pca_matrix_processing(actual_npcs = n_pcs, method = "irlba")
    cases[["canard_fixture_irlba_no_umap"]] <- list(
        matrix_processing = irlba_processor,
        compute_umap_embeddings = FALSE
    )
    if (env_flag("CA_CANARD_INCLUDE_UMAP", TRUE)) {
        cases[["canard_fixture_irlba_with_umap"]] <- list(
            matrix_processing = irlba_processor,
            compute_umap_embeddings = TRUE
        )
    }
}

results <- lapply(names(cases), function(label) {
    case_args <- c(common_arguments, cases[[label]])
    run_case(label, {
        suppressWarnings(do.call(automatic_stability_assessment, case_args))
    }, profile_dir = profile_dir, profile_interval = profile_interval)
})
results <- do.call(rbind, results)

print(results, row.names = FALSE)
cat(
    "fixture_app_dir=", app_dir,
    " fixture_genes=", nrow(expression_matrix),
    " fixture_cells=", ncol(expression_matrix),
    " steps=", paste(valid_steps, collapse = ","),
    " n_neigh=", paste(n_neigh_sequence, collapse = ","),
    " resolution=", paste(resolution_sequence, collapse = ","),
    " repetitions=", n_repetitions,
    "\n",
    sep = ""
)

output <- Sys.getenv("CA_CANARD_OUTPUT", "")
if (nzchar(output)) {
    utils::write.csv(results, output, row.names = FALSE)
}
